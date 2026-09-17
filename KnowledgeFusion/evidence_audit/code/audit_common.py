"""Neutral request, validation, and conversion utilities for EvidenceAudit.

The module is deliberately client-agnostic.  It prepares JSONL requests for a
structured-output LLM auditor, validates returned JSONL records, and converts
accepted records into the frozen CSV interface consumed by KnowledgeFusion.
It never calls a model or a network service.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PROMPT_DIR = PACKAGE_ROOT / "prompts"
SCHEMA_VERSION = "1.0.0"
PROTOCOL_VERSION = "7.0.0"

PROMPT_FILES = {
    "atomization": "01_context_aware_atomization.md",
    "mapping": "02_same_run_claim_evidence_mapping.md",
    "classification": "03_controlled_evidence_classification.md",
    "recheck": "04_revised_claim_recheck.md",
}

DESIGN_DIMENSIONS = {
    "uor_activity",
    "pba_structural_compatibility",
    "wet_chemical_synthesis",
    "alkaline_stability_corrosion",
    "substitution_and_synergy",
    "cost_environment",
}
MAPPING_STATUSES = {"related", "possibly_related", "unrelated"}
RELEVANCE = {"yes", "no"}
DIRECTIONS = {"support", "oppose", "mixed", "not_addressed"}
ENTAILMENTS = {"full", "partial", "contextual", "overreach", "none"}
SYSTEM_RELEVANCE = {
    "direct_PBA_UOR",
    "PBA_only",
    "UOR_only",
    "substitution",
    "foundational",
}
STUDY_TYPES = {"experimental", "computational", "review", "general_chemistry"}
TRACEABILITY = {
    "exact_chunk_DOI_UT",
    "abstract_DOI_UT",
    "identified_record",
    "partition_only",
    "untraceable",
}

CLAIM_FIELDS = [
    "governing_heading_excerpt",
    "context_excerpt",
    "claim_excerpt",
    "self_contained_claim",
    "subject",
    "object",
    "element_or_relation",
    "functional_role",
    "design_dimension",
    "claim_scope",
    "claim_status",
    "atomization_reason",
]
CLASSIFICATION_FIELDS = [
    "mapping_id",
    "relevant",
    "direction",
    "entailment",
    "system_relevance",
    "study_type",
    "traceability",
    "verified_excerpt",
    "classification_reason",
]
FROZEN_CLAIM_FIELDS = [
    "record_id",
    "candidate",
    "design_dimension",
    "claim_id",
    "included_after_human_review",
]
FROZEN_EVIDENCE_FIELDS = [
    "record_id",
    "candidate",
    "design_dimension",
    "independent_source_id",
    "source_packet_id",
    "DOI",
    "UT",
    "publication_year",
    "system_relevance",
    "study_type",
    "traceability",
    "entailment",
    "direction",
    "semantic_review_multiplier",
    "semantic_review_status",
    "included_in_ds",
]


class AuditError(ValueError):
    """Raised when an audit packet violates the public protocol contract."""


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_id(prefix: str, *parts: Any) -> str:
    payload = "\x1f".join(str(part) for part in parts)
    return f"{prefix}_{sha256_text(payload)[:16]}"


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            try:
                value = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise AuditError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            if not isinstance(value, dict):
                raise AuditError(f"{path}:{line_number}: each JSONL row must be an object")
            rows.append(value)
    return rows


def write_jsonl(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(canonical_json(dict(row)) + "\n")


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def require_fields(value: Mapping[str, Any], fields: Sequence[str], context: str) -> None:
    missing = [field for field in fields if field not in value]
    if missing:
        raise AuditError(f"{context}: missing fields: {', '.join(missing)}")


def require_string(value: Any, context: str, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise AuditError(f"{context}: expected a string")
    if not allow_empty and not value.strip():
        raise AuditError(f"{context}: value must not be empty")
    return value


def require_enum(value: Any, allowed: set[str], context: str) -> str:
    text = require_string(value, context)
    if text not in allowed:
        raise AuditError(f"{context}: {text!r} is not in {sorted(allowed)}")
    return text


def normalized_newlines(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def require_excerpt(excerpt: Any, full_text: str, context: str, allow_empty: bool = False) -> str:
    text = require_string(excerpt, context, allow_empty=allow_empty)
    if not text and allow_empty:
        return text
    if normalized_newlines(text) not in normalized_newlines(full_text):
        raise AuditError(f"{context}: excerpt is not a continuous span of the supplied text")
    return text


def prompt_metadata(task: str) -> dict[str, str]:
    if task not in PROMPT_FILES:
        raise AuditError(f"unknown audit task: {task}")
    path = PROMPT_DIR / PROMPT_FILES[task]
    text = path.read_text(encoding="utf-8")
    version_match = re.search(r"^Prompt-Version:[ \t]*(\S+)[ \t]*$", text, re.MULTILINE)
    id_match = re.search(r"^Prompt-ID:[ \t]*(\S+)[ \t]*$", text, re.MULTILINE)
    if not version_match or not id_match:
        raise AuditError(f"prompt metadata is incomplete: {path}")
    return {
        "prompt_id": id_match.group(1),
        "prompt_version": version_match.group(1),
        "prompt_sha256": sha256_file(path),
    }


def request_digest(request: Mapping[str, Any]) -> str:
    payload = {key: value for key, value in request.items() if key != "request_hash"}
    return sha256_text(canonical_json(payload))


def make_request(task: str, input_payload: Mapping[str, Any]) -> dict[str, Any]:
    meta = prompt_metadata(task)
    identity = sha256_text(canonical_json(input_payload))[:16]
    request = {
        "schema_version": SCHEMA_VERSION,
        "task": task,
        "protocol_version": PROTOCOL_VERSION,
        "request_id": f"request_{task}_{identity}",
        **meta,
        "input": dict(input_payload),
    }
    request["request_hash"] = request_digest(request)
    return request


def validate_request_record(request: Mapping[str, Any]) -> None:
    require_fields(
        request,
        [
            "schema_version",
            "task",
            "protocol_version",
            "request_id",
            "request_hash",
            "prompt_id",
            "prompt_version",
            "prompt_sha256",
            "input",
        ],
        "request",
    )
    task = require_enum(request["task"], set(PROMPT_FILES), "request.task")
    expected_meta = prompt_metadata(task)
    for field, expected in expected_meta.items():
        if request[field] != expected:
            raise AuditError(f"request.{field}: expected {expected!r}")
    if request["protocol_version"] != PROTOCOL_VERSION:
        raise AuditError("request.protocol_version does not match the public protocol")
    if request["request_hash"] != request_digest(request):
        raise AuditError("request.request_hash does not match the request content")
    if not isinstance(request["input"], dict):
        raise AuditError("request.input must be an object")


def read_run_packets(path: Path) -> list[dict[str, Any]]:
    packets = load_jsonl(path)
    seen_runs: set[str] = set()
    for index, packet in enumerate(packets):
        context = f"run packet {index + 1}"
        require_fields(packet, ["record_id", "candidate", "final_answer", "source_packets"], context)
        record_id = require_string(packet["record_id"], f"{context}.record_id")
        require_string(packet["candidate"], f"{context}.candidate")
        require_string(packet["final_answer"], f"{context}.final_answer")
        if record_id in seen_runs:
            raise AuditError(f"duplicate record_id: {record_id}")
        seen_runs.add(record_id)
        if not isinstance(packet["source_packets"], list):
            raise AuditError(f"{context}.source_packets must be an array")
        source_ids: set[str] = set()
        for source_index, source in enumerate(packet["source_packets"]):
            source_context = f"{context}.source_packets[{source_index}]"
            require_fields(source, ["source_packet_id", "evidence_text"], source_context)
            source_id = require_string(source["source_packet_id"], f"{source_context}.source_packet_id")
            require_string(source["evidence_text"], f"{source_context}.evidence_text")
            if source_id in source_ids:
                raise AuditError(f"{context}: duplicate source_packet_id: {source_id}")
            source_ids.add(source_id)
    return packets


def index_run_packets(packets: Sequence[Mapping[str, Any]]) -> dict[str, Mapping[str, Any]]:
    return {str(packet["record_id"]): packet for packet in packets}


def source_index(packet: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    return {str(source["source_packet_id"]): source for source in packet["source_packets"]}


def prepare_atomization_requests(run_packets: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    requests: list[dict[str, Any]] = []
    for packet in run_packets:
        requests.append(
            make_request(
                "atomization",
                {
                    "record_id": packet["record_id"],
                    "candidate": packet["candidate"],
                    "final_answer": packet["final_answer"],
                    "allowed_functional_roles": packet.get("allowed_functional_roles", []),
                    "allowed_design_dimensions": packet.get(
                        "allowed_design_dimensions", sorted(DESIGN_DIMENSIONS)
                    ),
                },
            )
        )
    return requests


def _validated_claims_by_run(rows: Sequence[Mapping[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}
    for response in rows:
        for claim in response.get("claims", []):
            result.setdefault(str(claim["record_id"]), []).append(dict(claim))
    return result


def prepare_mapping_requests(
    run_packets: Sequence[Mapping[str, Any]],
    atomization_responses: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    claims_by_run = _validated_claims_by_run(atomization_responses)
    requests: list[dict[str, Any]] = []
    for packet in run_packets:
        record_id = str(packet["record_id"])
        claims = claims_by_run.get(record_id, [])
        if not claims:
            raise AuditError(f"no validated atomized claims found for {record_id}")
        compact_claims = [
            {
                "claim_id": claim["claim_id"],
                "self_contained_claim": claim["self_contained_claim"],
                "design_dimension": claim["design_dimension"],
            }
            for claim in claims
        ]
        for source in packet["source_packets"]:
            requests.append(
                make_request(
                    "mapping",
                    {
                        "record_id": record_id,
                        "candidate": packet["candidate"],
                        "source_packet_id": source["source_packet_id"],
                        "source_metadata": {
                            "DOI": source.get("DOI", ""),
                            "UT": source.get("UT", ""),
                            "publication_year": source.get(
                                "publication_year", source.get("study_year", "")
                            ),
                            "independent_source_id": source.get(
                                "independent_source_id", source.get("independence_group", "")
                            ),
                        },
                        "evidence_text": source["evidence_text"],
                        "claims": compact_claims,
                    },
                )
            )
    return requests


def prepare_classification_requests(
    run_packets: Sequence[Mapping[str, Any]],
    mapping_responses: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    runs = index_run_packets(run_packets)
    requests: list[dict[str, Any]] = []
    for response in mapping_responses:
        for mapping in response.get("mappings", []):
            if mapping["mapping_status"] == "unrelated":
                continue
            record_id = str(mapping["record_id"])
            if record_id not in runs:
                raise AuditError(f"mapping references unknown record_id: {record_id}")
            packet = runs[record_id]
            sources = source_index(packet)
            source_id = str(mapping["source_packet_id"])
            if source_id not in sources:
                raise AuditError(f"mapping references a source outside its run: {source_id}")
            source = sources[source_id]
            requests.append(
                make_request(
                    "classification",
                    {
                        "record_id": record_id,
                        "candidate": packet["candidate"],
                        "claim_id": mapping["claim_id"],
                        "self_contained_claim": mapping["self_contained_claim"],
                        "design_dimension": mapping["design_dimension"],
                        "mapping_id": mapping["mapping_id"],
                        "mapping_status": mapping["mapping_status"],
                        "mapping_excerpt": mapping["mapping_excerpt"],
                        "source_packet_id": source_id,
                        "source_metadata": {
                            "DOI": source.get("DOI", ""),
                            "UT": source.get("UT", ""),
                            "publication_year": source.get(
                                "publication_year", source.get("study_year", "")
                            ),
                            "independent_source_id": source.get(
                                "independent_source_id", source.get("independence_group", "")
                            ),
                        },
                        "evidence_text": source["evidence_text"],
                    },
                )
            )
    return requests


def prepare_recheck_requests(
    run_packets: Sequence[Mapping[str, Any]], revised_claims: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    runs = index_run_packets(run_packets)
    requests: list[dict[str, Any]] = []
    for revised in revised_claims:
        require_fields(
            revised,
            ["record_id", "claim_id", "self_contained_claim", "design_dimension"],
            "revised claim",
        )
        record_id = str(revised["record_id"])
        if record_id not in runs:
            raise AuditError(f"revised claim references unknown record_id: {record_id}")
        packet = runs[record_id]
        requests.append(
            make_request(
                "recheck",
                {
                    "record_id": record_id,
                    "candidate": packet["candidate"],
                    "claim_id": revised["claim_id"],
                    "self_contained_claim": revised["self_contained_claim"],
                    "design_dimension": revised["design_dimension"],
                    "source_packets": packet["source_packets"],
                },
            )
        )
    return requests


def prepare_requests(
    phase: str,
    run_packets_path: Path,
    output_path: Path,
    responses_path: Path | None = None,
    revised_claims_path: Path | None = None,
) -> list[dict[str, Any]]:
    packets = read_run_packets(run_packets_path)
    if phase == "atomization":
        requests = prepare_atomization_requests(packets)
    elif phase == "mapping":
        if responses_path is None:
            raise AuditError("mapping preparation requires --responses with validated atomization JSONL")
        requests = prepare_mapping_requests(packets, load_jsonl(responses_path))
    elif phase == "classification":
        if responses_path is None:
            raise AuditError("classification preparation requires --responses with validated mapping JSONL")
        requests = prepare_classification_requests(packets, load_jsonl(responses_path))
    elif phase == "recheck":
        if revised_claims_path is None:
            raise AuditError("recheck preparation requires --revised-claims")
        requests = prepare_recheck_requests(packets, load_jsonl(revised_claims_path))
    else:
        raise AuditError(f"unsupported phase: {phase}")
    write_jsonl(output_path, requests)
    return requests


def _validate_response_metadata(response: Mapping[str, Any], request: Mapping[str, Any]) -> None:
    required = [
        "request_id",
        "request_hash",
        "prompt_id",
        "prompt_sha256",
        "protocol_version",
        "review_model",
        "review_date",
    ]
    require_fields(response, required, f"response {request['request_id']}")
    for field in [
        "request_id",
        "request_hash",
        "prompt_id",
        "prompt_sha256",
        "protocol_version",
    ]:
        if response[field] != request[field]:
            raise AuditError(f"response {request['request_id']}: {field} does not match its request")
    require_string(response["review_model"], "response.review_model")
    date = require_string(response["review_date"], "response.review_date")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date):
        raise AuditError("response.review_date must use YYYY-MM-DD")


def _validate_atomization(response: Mapping[str, Any], request: Mapping[str, Any]) -> dict[str, Any]:
    claims = response.get("claims")
    if not isinstance(claims, list):
        raise AuditError("atomization response.claims must be an array")
    final_answer = str(request["input"]["final_answer"])
    result = dict(response)
    validated: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, raw_claim in enumerate(claims):
        if not isinstance(raw_claim, dict):
            raise AuditError(f"claims[{index}] must be an object")
        require_fields(raw_claim, CLAIM_FIELDS, f"claims[{index}]")
        claim = dict(raw_claim)
        context_excerpt = require_excerpt(
            claim["context_excerpt"], final_answer, f"claims[{index}].context_excerpt"
        )
        require_excerpt(claim["claim_excerpt"], context_excerpt, f"claims[{index}].claim_excerpt")
        heading = require_string(
            claim["governing_heading_excerpt"], f"claims[{index}].governing_heading_excerpt"
        )
        if heading != "NO_EXPLICIT_HEADING":
            require_excerpt(heading, final_answer, f"claims[{index}].governing_heading_excerpt")
        for field in [
            "self_contained_claim",
            "subject",
            "object",
            "element_or_relation",
            "functional_role",
            "claim_scope",
            "claim_status",
            "atomization_reason",
        ]:
            require_string(claim[field], f"claims[{index}].{field}")
        require_enum(claim["design_dimension"], DESIGN_DIMENSIONS, f"claims[{index}].design_dimension")
        claim_id = stable_id(
            "claim",
            request["input"]["record_id"],
            index,
            claim["self_contained_claim"],
            claim["claim_excerpt"],
        )
        if claim_id in seen_ids:
            raise AuditError(f"duplicate deterministic claim_id: {claim_id}")
        seen_ids.add(claim_id)
        claim.update(
            {
                "claim_id": claim_id,
                "record_id": request["input"]["record_id"],
                "candidate": request["input"]["candidate"],
            }
        )
        validated.append(claim)
    result["claims"] = validated
    return result


def _validate_mapping(response: Mapping[str, Any], request: Mapping[str, Any]) -> dict[str, Any]:
    mappings = response.get("mappings")
    if not isinstance(mappings, list):
        raise AuditError("mapping response.mappings must be an array")
    input_claims = {str(claim["claim_id"]): claim for claim in request["input"]["claims"]}
    if len(mappings) != len(input_claims):
        raise AuditError("mapping response must contain exactly one item per supplied claim_id")
    seen: set[str] = set()
    validated: list[dict[str, Any]] = []
    evidence_text = str(request["input"]["evidence_text"])
    for index, raw_mapping in enumerate(mappings):
        if not isinstance(raw_mapping, dict):
            raise AuditError(f"mappings[{index}] must be an object")
        require_fields(
            raw_mapping,
            ["claim_id", "mapping_status", "mapping_excerpt", "mapping_reason"],
            f"mappings[{index}]",
        )
        claim_id = require_string(raw_mapping["claim_id"], f"mappings[{index}].claim_id")
        if claim_id not in input_claims or claim_id in seen:
            raise AuditError(f"mappings[{index}]: unknown or duplicate claim_id {claim_id}")
        seen.add(claim_id)
        status = require_enum(
            raw_mapping["mapping_status"], MAPPING_STATUSES, f"mappings[{index}].mapping_status"
        )
        excerpt = require_string(
            raw_mapping["mapping_excerpt"], f"mappings[{index}].mapping_excerpt", allow_empty=True
        )
        if status == "unrelated":
            if excerpt:
                raise AuditError("unrelated mappings must have an empty mapping_excerpt")
        else:
            require_excerpt(excerpt, evidence_text, f"mappings[{index}].mapping_excerpt")
        require_string(raw_mapping["mapping_reason"], f"mappings[{index}].mapping_reason")
        claim = input_claims[claim_id]
        mapping = dict(raw_mapping)
        mapping.update(
            {
                "mapping_id": stable_id(
                    "mapping",
                    request["input"]["record_id"],
                    request["input"]["source_packet_id"],
                    claim_id,
                ),
                "record_id": request["input"]["record_id"],
                "candidate": request["input"]["candidate"],
                "source_packet_id": request["input"]["source_packet_id"],
                "source_metadata": request["input"].get("source_metadata", {}),
                "self_contained_claim": claim["self_contained_claim"],
                "design_dimension": claim["design_dimension"],
            }
        )
        validated.append(mapping)
    result = dict(response)
    result["mappings"] = validated
    return result


def _validate_classification_item(
    item: Mapping[str, Any], evidence_text: str, expected_mapping_id: str, context: str
) -> dict[str, Any]:
    require_fields(item, CLASSIFICATION_FIELDS, context)
    if item["mapping_id"] != expected_mapping_id:
        raise AuditError(f"{context}.mapping_id does not match the request")
    relevant = require_enum(item["relevant"], RELEVANCE, f"{context}.relevant")
    direction = require_enum(item["direction"], DIRECTIONS, f"{context}.direction")
    entailment = require_enum(item["entailment"], ENTAILMENTS, f"{context}.entailment")
    require_enum(item["system_relevance"], SYSTEM_RELEVANCE, f"{context}.system_relevance")
    require_enum(item["study_type"], STUDY_TYPES, f"{context}.study_type")
    require_enum(item["traceability"], TRACEABILITY, f"{context}.traceability")
    verified = require_string(item["verified_excerpt"], f"{context}.verified_excerpt", allow_empty=True)
    require_string(item["classification_reason"], f"{context}.classification_reason")
    if direction == "not_addressed" and entailment != "none":
        raise AuditError(f"{context}: not_addressed requires entailment=none")
    if relevant == "no" and direction != "not_addressed":
        raise AuditError(f"{context}: relevant=no requires direction=not_addressed")
    if direction in {"support", "oppose", "mixed"}:
        require_excerpt(verified, evidence_text, f"{context}.verified_excerpt")
    elif verified:
        require_excerpt(verified, evidence_text, f"{context}.verified_excerpt")
    return dict(item)


def _validate_classification(
    response: Mapping[str, Any], request: Mapping[str, Any]
) -> dict[str, Any]:
    classifications = response.get("classifications")
    if not isinstance(classifications, list) or len(classifications) != 1:
        raise AuditError("classification response must contain exactly one classification")
    item = _validate_classification_item(
        classifications[0],
        str(request["input"]["evidence_text"]),
        str(request["input"]["mapping_id"]),
        "classifications[0]",
    )
    item.update(
        {
            "record_id": request["input"]["record_id"],
            "candidate": request["input"]["candidate"],
            "claim_id": request["input"]["claim_id"],
            "self_contained_claim": request["input"]["self_contained_claim"],
            "design_dimension": request["input"]["design_dimension"],
            "source_packet_id": request["input"]["source_packet_id"],
            "source_metadata": request["input"].get("source_metadata", {}),
        }
    )
    result = dict(response)
    result["classifications"] = [item]
    return result


def _validate_recheck(response: Mapping[str, Any], request: Mapping[str, Any]) -> dict[str, Any]:
    assessments = response.get("source_assessments")
    if not isinstance(assessments, list):
        raise AuditError("recheck response.source_assessments must be an array")
    sources = {str(item["source_packet_id"]): item for item in request["input"]["source_packets"]}
    if len(assessments) != len(sources):
        raise AuditError("recheck response must contain exactly one assessment per source_packet_id")
    seen: set[str] = set()
    validated: list[dict[str, Any]] = []
    for index, raw in enumerate(assessments):
        context = f"source_assessments[{index}]"
        require_fields(
            raw,
            ["source_packet_id", "mapping_status", "mapping_excerpt", "mapping_reason"],
            context,
        )
        source_id = require_string(raw["source_packet_id"], f"{context}.source_packet_id")
        if source_id not in sources or source_id in seen:
            raise AuditError(f"{context}: unknown or duplicate source_packet_id {source_id}")
        seen.add(source_id)
        source = sources[source_id]
        status = require_enum(raw["mapping_status"], MAPPING_STATUSES, f"{context}.mapping_status")
        excerpt = require_string(raw["mapping_excerpt"], f"{context}.mapping_excerpt", allow_empty=True)
        if status == "unrelated":
            if excerpt:
                raise AuditError(f"{context}: unrelated source must have empty mapping_excerpt")
            for field in CLASSIFICATION_FIELDS[1:]:
                if raw.get(field, "") not in ("", None):
                    raise AuditError(f"{context}: unrelated source must leave {field} empty")
        else:
            require_excerpt(excerpt, str(source["evidence_text"]), f"{context}.mapping_excerpt")
            classification = dict(raw)
            classification["mapping_id"] = stable_id(
                "mapping", request["input"]["record_id"], source_id, request["input"]["claim_id"]
            )
            _validate_classification_item(
                classification,
                str(source["evidence_text"]),
                classification["mapping_id"],
                context,
            )
        item = dict(raw)
        item["mapping_id"] = stable_id(
            "mapping", request["input"]["record_id"], source_id, request["input"]["claim_id"]
        )
        item.update(
            {
                "record_id": request["input"]["record_id"],
                "candidate": request["input"]["candidate"],
                "claim_id": request["input"]["claim_id"],
                "self_contained_claim": request["input"]["self_contained_claim"],
                "design_dimension": request["input"]["design_dimension"],
                "source_metadata": {
                    "DOI": source.get("DOI", ""),
                    "UT": source.get("UT", ""),
                    "publication_year": source.get(
                        "publication_year", source.get("study_year", "")
                    ),
                    "independent_source_id": source.get(
                        "independent_source_id", source.get("independence_group", "")
                    ),
                },
            }
        )
        validated.append(item)
    result = dict(response)
    result["source_assessments"] = validated
    return result


def validate_responses(
    requests_path: Path, responses_path: Path, output_path: Path
) -> list[dict[str, Any]]:
    requests = load_jsonl(requests_path)
    responses = load_jsonl(responses_path)
    request_by_id: dict[str, dict[str, Any]] = {}
    for request in requests:
        validate_request_record(request)
        request_id = str(request["request_id"])
        if request_id in request_by_id:
            raise AuditError(f"duplicate request_id: {request_id}")
        request_by_id[request_id] = request
    response_by_id: dict[str, dict[str, Any]] = {}
    for response in responses:
        request_id = require_string(response.get("request_id"), "response.request_id")
        if request_id not in request_by_id:
            raise AuditError(f"response references unknown request_id: {request_id}")
        if request_id in response_by_id:
            raise AuditError(f"duplicate response for request_id: {request_id}")
        response_by_id[request_id] = response
    missing = sorted(set(request_by_id) - set(response_by_id))
    if missing:
        raise AuditError(f"missing responses for {len(missing)} request(s): {', '.join(missing[:5])}")
    validated: list[dict[str, Any]] = []
    for request in requests:
        response = response_by_id[str(request["request_id"])]
        _validate_response_metadata(response, request)
        task = str(request["task"])
        if task == "atomization":
            item = _validate_atomization(response, request)
        elif task == "mapping":
            item = _validate_mapping(response, request)
        elif task == "classification":
            item = _validate_classification(response, request)
        elif task == "recheck":
            item = _validate_recheck(response, request)
        else:
            raise AuditError(f"unsupported request task: {task}")
        item["validated_request_id"] = request["request_id"]
        validated.append(item)
    write_jsonl(output_path, validated)
    return validated


def independent_source_id(metadata: Mapping[str, Any], source_packet_id: str) -> str:
    declared = str(metadata.get("independent_source_id", "")).strip()
    if declared:
        return declared
    doi = str(metadata.get("DOI", "")).strip().lower()
    if doi:
        return f"doi:{doi}"
    ut = str(metadata.get("UT", "")).strip()
    if ut:
        return f"ut:{ut}"
    return f"packet:{source_packet_id}"


def _flatten_classifications(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for row in rows:
        for item in row.get("classifications", []):
            result.append(dict(item))
        for item in row.get("source_assessments", []):
            if item.get("mapping_status") != "unrelated":
                result.append(dict(item))
    return result


def _flatten_mappings(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for row in rows:
        for item in row.get("mappings", []):
            if item.get("mapping_status") != "unrelated":
                result.append(dict(item))
        for item in row.get("source_assessments", []):
            if item.get("mapping_status") != "unrelated":
                result.append(dict(item))
    return result


def build_frozen_inputs(
    atomization_path: Path,
    mapping_path: Path,
    classification_path: Path,
    adjudication_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    atomization = load_jsonl(atomization_path)
    mapping = load_jsonl(mapping_path)
    classification = load_jsonl(classification_path)
    adjudication = load_json(adjudication_path)
    require_fields(
        adjudication,
        ["schema_version", "semantic_discount", "claim_decisions", "source_packet_decisions"],
        "adjudication",
    )
    discount = float(adjudication["semantic_discount"])
    if not 0.0 <= discount <= 1.0:
        raise AuditError("adjudication.semantic_discount must be between 0 and 1")

    claims = [claim for row in atomization for claim in row.get("claims", [])]
    claim_by_id = {str(claim["claim_id"]): claim for claim in claims}
    if len(claim_by_id) != len(claims):
        raise AuditError("atomization input contains duplicate claim_id values")
    claim_decisions: dict[str, Mapping[str, Any]] = {}
    for item in adjudication["claim_decisions"]:
        require_fields(item, ["claim_id", "decision"], "claim decision")
        claim_id = str(item["claim_id"])
        if claim_id not in claim_by_id or claim_id in claim_decisions:
            raise AuditError(f"unknown or duplicate claim decision: {claim_id}")
        require_enum(item["decision"], {"include", "exclude", "revise"}, "claim decision.decision")
        claim_decisions[claim_id] = item
    missing_claims = sorted(set(claim_by_id) - set(claim_decisions))
    if missing_claims:
        raise AuditError(f"adjudication is missing {len(missing_claims)} claim decision(s)")

    source_decisions: dict[tuple[str, str], Mapping[str, Any]] = {}
    for item in adjudication["source_packet_decisions"]:
        require_fields(item, ["record_id", "source_packet_id", "decision"], "source decision")
        key = (str(item["record_id"]), str(item["source_packet_id"]))
        if key in source_decisions:
            raise AuditError(f"duplicate source decision: {key}")
        require_enum(
            item["decision"],
            {"accepted", "machine_unreviewed", "revised", "rejected"},
            "source decision.decision",
        )
        source_decisions[key] = item

    frozen_claims: list[dict[str, Any]] = []
    included_claim_ids: set[str] = set()
    for claim_id, claim in claim_by_id.items():
        decision = claim_decisions[claim_id]
        if decision["decision"] == "exclude":
            continue
        design_dimension = str(decision.get("design_dimension", claim["design_dimension"]))
        require_enum(design_dimension, DESIGN_DIMENSIONS, "claim decision.design_dimension")
        included_claim_ids.add(claim_id)
        frozen_claims.append(
            {
                "record_id": claim["record_id"],
                "candidate": claim["candidate"],
                "design_dimension": design_dimension,
                "claim_id": claim_id,
                "included_after_human_review": True,
            }
        )

    mappings = _flatten_mappings(mapping)
    mapping_by_id: dict[str, dict[str, Any]] = {}
    for item in mappings:
        mapping_id = str(item["mapping_id"])
        if mapping_id in mapping_by_id:
            raise AuditError(f"mapping input contains duplicate mapping_id: {mapping_id}")
        mapping_by_id[mapping_id] = item

    classifications = _flatten_classifications(classification)
    classified_ids: set[str] = set()
    for item in classifications:
        mapping_id = str(item["mapping_id"])
        if mapping_id not in mapping_by_id or mapping_id in classified_ids:
            raise AuditError(f"classification has an unknown or duplicate mapping_id: {mapping_id}")
        classified_ids.add(mapping_id)
        mapped = mapping_by_id[mapping_id]
        for field in ["record_id", "claim_id", "source_packet_id"]:
            if str(item[field]) != str(mapped[field]):
                raise AuditError(f"classification {mapping_id} disagrees with its mapping on {field}")
    missing_classifications = sorted(set(mapping_by_id) - classified_ids)
    if missing_classifications:
        raise AuditError(
            f"classification input is missing {len(missing_classifications)} mapped pair(s)"
        )

    frozen_evidence: list[dict[str, Any]] = []
    for item in classifications:
        claim_id = str(item["claim_id"])
        if claim_id not in included_claim_ids:
            continue
        record_id = str(item["record_id"])
        source_packet_id = str(item["source_packet_id"])
        key = (record_id, source_packet_id)
        if key not in source_decisions:
            raise AuditError(f"adjudication is missing a source decision for {key}")
        source_decision = str(source_decisions[key]["decision"])
        if source_decision == "rejected":
            continue
        status = "machine_unreviewed" if source_decision == "machine_unreviewed" else "human_accepted"
        multiplier = discount if status == "machine_unreviewed" else 1.0
        metadata = item.get("source_metadata", {})
        frozen_evidence.append(
            {
                "record_id": record_id,
                "candidate": item["candidate"],
                "design_dimension": claim_by_id[claim_id]["design_dimension"],
                "independent_source_id": independent_source_id(metadata, source_packet_id),
                "source_packet_id": source_packet_id,
                "DOI": metadata.get("DOI", ""),
                "UT": metadata.get("UT", ""),
                "publication_year": metadata.get("publication_year", ""),
                "system_relevance": item["system_relevance"],
                "study_type": item["study_type"],
                "traceability": item["traceability"],
                "entailment": item["entailment"],
                "direction": item["direction"],
                "semantic_review_multiplier": format(multiplier, ".16g"),
                "semantic_review_status": status,
                "included_in_ds": "yes",
            }
        )

    frozen_claims.sort(key=lambda row: (row["record_id"], row["claim_id"]))
    frozen_evidence.sort(
        key=lambda row: (
            row["record_id"],
            row["source_packet_id"],
            row["design_dimension"],
            row["independent_source_id"],
        )
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "frozen_claims.csv", FROZEN_CLAIM_FIELDS, frozen_claims)
    write_csv(output_dir / "frozen_evidence.csv", FROZEN_EVIDENCE_FIELDS, frozen_evidence)

    source_counts = Counter(str(value["decision"]) for value in source_decisions.values())
    summary = {
        "accepted_auditor_packets": source_counts["accepted"] + source_counts["revised"],
        "auditor_only_source_packets": int(
            adjudication.get(
                "auditor_only_source_packets",
                source_counts["accepted"] + source_counts["revised"],
            )
        ),
        "frozen_claim_count": len(frozen_claims),
        "frozen_evidence_count": len(frozen_evidence),
        "pending_actions": 0,
        "ready_for_fusion": True,
        "rejected_packets": source_counts["rejected"],
        "review_coverage": adjudication.get("review_coverage", []),
        "revised_packets": source_counts["revised"],
        "schema_version": SCHEMA_VERSION,
        "semantic_discount_for_unreviewed_sources": discount,
    }
    write_json(output_dir / "frozen_review_summary.json", summary)
    return summary


def verify_reference_frozen_inputs(input_dir: Path) -> dict[str, Any]:
    claims = load_csv(input_dir / "frozen_claims.csv")
    evidence = load_csv(input_dir / "frozen_evidence.csv")
    summary = load_json(input_dir / "frozen_review_summary.json")
    if list(claims[0]) != FROZEN_CLAIM_FIELDS:
        raise AuditError("reference frozen_claims.csv fields do not match the public interface")
    if list(evidence[0]) != FROZEN_EVIDENCE_FIELDS:
        raise AuditError("reference frozen_evidence.csv fields do not match the public interface")
    if len(claims) != int(summary["frozen_claim_count"]):
        raise AuditError("reference frozen claim count disagrees with its summary")
    if len(evidence) != int(summary["frozen_evidence_count"]):
        raise AuditError("reference frozen evidence count disagrees with its summary")
    multipliers = {
        row["semantic_review_multiplier"]
        for row in evidence
        if row["semantic_review_status"] == "machine_unreviewed"
    }
    expected_discount = str(summary["semantic_discount_for_unreviewed_sources"])
    if multipliers != {expected_discount}:
        raise AuditError("reference unreviewed evidence multiplier disagrees with its summary")
    return {
        "frozen_claim_count": len(claims),
        "frozen_evidence_count": len(evidence),
        "accepted_auditor_packets": int(summary["accepted_auditor_packets"]),
        "semantic_discount_for_unreviewed_sources": float(expected_discount),
    }


def path_arg(value: str) -> Path:
    return Path(value).expanduser().resolve()


def add_common_validation_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--requests", type=path_arg, required=True)
    parser.add_argument("--responses", type=path_arg, required=True)
    parser.add_argument("--output", type=path_arg, required=True)
