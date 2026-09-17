"""End-to-end tests for the neutral EvidenceAudit interface."""

from __future__ import annotations

import argparse
import copy
import tempfile
import unittest
from pathlib import Path

from audit_common import (
    PACKAGE_ROOT,
    AuditError,
    build_frozen_inputs,
    prepare_atomization_requests,
    prepare_classification_requests,
    prepare_mapping_requests,
    prepare_recheck_requests,
    sha256_file,
    validate_responses,
    verify_reference_frozen_inputs,
    write_json,
    write_jsonl,
)


EXPECTED_PROMPT_HASHES = {
    "01_context_aware_atomization.md": "f5d92fe97bfd1c41c61dca974e6ed3bc8041297d8d92854faa43942d10018dae",
    "02_same_run_claim_evidence_mapping.md": "8bfef05f5f838b25c65abc2f0eb5cd17a39b3b0a07c62d646444cbcfd4d1b184",
    "03_controlled_evidence_classification.md": "2610f029293609b9125ba700597d22be01ae73b627d2e642cb5e485c2787c8dc",
    "04_revised_claim_recheck.md": "06f79af345d58ae76a281eb11569bcc3538e17b4864e5212e2a3013b6733944d",
}


def response_header(request: dict) -> dict:
    return {
        "request_id": request["request_id"],
        "request_hash": request["request_hash"],
        "prompt_id": request["prompt_id"],
        "prompt_sha256": request["prompt_sha256"],
        "protocol_version": request["protocol_version"],
        "review_model": "generic-structured-output-model",
        "review_date": "2000-01-01",
    }


def demo_packet() -> dict:
    return {
        "record_id": "run_demo_001",
        "candidate": "Fe-Co-Ni-Mn-Zn",
        "final_answer": (
            "## Nickel role\n"
            "Nickel may provide redox-active sites under alkaline UOR conditions."
        ),
        "allowed_functional_roles": ["redox_activity"],
        "allowed_design_dimensions": ["uor_activity"],
        "source_packets": [
            {
                "source_packet_id": "source_demo_001",
                "evidence_text": (
                    "Synthetic example source: nickel redox centers were examined "
                    "under alkaline UOR conditions."
                ),
                "DOI": "",
                "UT": "",
                "publication_year": "",
            }
        ],
    }


def run_demo(base: Path) -> dict:
    input_dir = base / "input"
    request_dir = base / "requests"
    response_dir = base / "responses"
    expected_dir = base / "expected"
    for directory in (input_dir, request_dir, response_dir, expected_dir):
        directory.mkdir(parents=True, exist_ok=True)

    packets = [demo_packet()]
    run_packet_path = input_dir / "run_packets.jsonl"
    write_jsonl(run_packet_path, packets)

    atom_requests = prepare_atomization_requests(packets)
    atom_request_path = request_dir / "atomization_requests.jsonl"
    write_jsonl(atom_request_path, atom_requests)
    atom_response = {
        **response_header(atom_requests[0]),
        "claims": [
            {
                "governing_heading_excerpt": "## Nickel role",
                "context_excerpt": (
                    "Nickel may provide redox-active sites under alkaline UOR conditions."
                ),
                "claim_excerpt": (
                    "Nickel may provide redox-active sites under alkaline UOR conditions."
                ),
                "self_contained_claim": (
                    "Nickel may provide redox-active sites under alkaline UOR conditions."
                ),
                "subject": "Nickel",
                "object": "redox-active sites",
                "element_or_relation": "Ni",
                "functional_role": "redox_activity",
                "design_dimension": "uor_activity",
                "claim_scope": "candidate component",
                "claim_status": "explicit_claim",
                "atomization_reason": "One explicit material-selection proposition.",
            }
        ],
    }
    atom_response_path = response_dir / "atomization_responses.jsonl"
    atom_validated_path = expected_dir / "validated_atomization.jsonl"
    write_jsonl(atom_response_path, [atom_response])
    atom_validated = validate_responses(
        atom_request_path, atom_response_path, atom_validated_path
    )

    mapping_requests = prepare_mapping_requests(packets, atom_validated)
    mapping_request_path = request_dir / "mapping_requests.jsonl"
    write_jsonl(mapping_request_path, mapping_requests)
    mapping_response = {
        **response_header(mapping_requests[0]),
        "mappings": [
            {
                "claim_id": mapping_requests[0]["input"]["claims"][0]["claim_id"],
                "mapping_status": "related",
                "mapping_excerpt": (
                    "nickel redox centers were examined under alkaline UOR conditions"
                ),
                "mapping_reason": "The source addresses the same element and UOR context.",
            }
        ],
    }
    mapping_response_path = response_dir / "mapping_responses.jsonl"
    mapping_validated_path = expected_dir / "validated_mapping.jsonl"
    write_jsonl(mapping_response_path, [mapping_response])
    mapping_validated = validate_responses(
        mapping_request_path, mapping_response_path, mapping_validated_path
    )

    classification_requests = prepare_classification_requests(packets, mapping_validated)
    classification_request_path = request_dir / "classification_requests.jsonl"
    write_jsonl(classification_request_path, classification_requests)
    classification_response = {
        **response_header(classification_requests[0]),
        "classifications": [
            {
                "mapping_id": classification_requests[0]["input"]["mapping_id"],
                "relevant": "yes",
                "direction": "support",
                "entailment": "partial",
                "system_relevance": "direct_PBA_UOR",
                "study_type": "experimental",
                "traceability": "partition_only",
                "verified_excerpt": (
                    "nickel redox centers were examined under alkaline UOR conditions"
                ),
                "classification_reason": (
                    "The synthetic source directly addresses the claimed redox role."
                ),
            }
        ],
    }
    classification_response_path = response_dir / "classification_responses.jsonl"
    classification_validated_path = expected_dir / "validated_classification.jsonl"
    write_jsonl(classification_response_path, [classification_response])
    validate_responses(
        classification_request_path,
        classification_response_path,
        classification_validated_path,
    )

    claim_id = atom_validated[0]["claims"][0]["claim_id"]
    revised_claims = [
        {
            "record_id": "run_demo_001",
            "claim_id": claim_id,
            "self_contained_claim": (
                "Nickel may contribute redox-active sites under alkaline UOR conditions."
            ),
            "design_dimension": "uor_activity",
        }
    ]
    revised_path = input_dir / "revised_claims.jsonl"
    write_jsonl(revised_path, revised_claims)
    recheck_requests = prepare_recheck_requests(packets, revised_claims)
    recheck_request_path = request_dir / "recheck_requests.jsonl"
    write_jsonl(recheck_request_path, recheck_requests)
    recheck_response = {
        **response_header(recheck_requests[0]),
        "source_assessments": [
            {
                "source_packet_id": "source_demo_001",
                "mapping_status": "related",
                "mapping_excerpt": (
                    "nickel redox centers were examined under alkaline UOR conditions"
                ),
                "mapping_reason": "The source addresses the revised claim.",
                "relevant": "yes",
                "direction": "support",
                "entailment": "partial",
                "system_relevance": "direct_PBA_UOR",
                "study_type": "experimental",
                "traceability": "partition_only",
                "verified_excerpt": (
                    "nickel redox centers were examined under alkaline UOR conditions"
                ),
                "classification_reason": "The source supports the narrower claim.",
            }
        ],
    }
    recheck_response_path = response_dir / "recheck_responses.jsonl"
    recheck_validated_path = expected_dir / "validated_recheck.jsonl"
    write_jsonl(recheck_response_path, [recheck_response])
    validate_responses(recheck_request_path, recheck_response_path, recheck_validated_path)

    adjudication = {
        "schema_version": "1.0.0",
        "semantic_discount": 0.9674014658411585,
        "auditor_only_source_packets": 1,
        "claim_decisions": [{"claim_id": claim_id, "decision": "include"}],
        "source_packet_decisions": [
            {
                "record_id": "run_demo_001",
                "source_packet_id": "source_demo_001",
                "decision": "accepted",
            }
        ],
        "review_coverage": [
            {"category": "synthetic_example", "population": 1, "reviewed": 1}
        ],
    }
    adjudication_path = input_dir / "adjudication.json"
    write_json(adjudication_path, adjudication)
    summary = build_frozen_inputs(
        atom_validated_path,
        mapping_validated_path,
        classification_validated_path,
        adjudication_path,
        expected_dir,
    )
    return {
        "summary": summary,
        "classification_request": classification_requests[0],
        "classification_response": classification_response,
        "classification_request_path": classification_request_path,
    }


class EvidenceAuditTests(unittest.TestCase):
    def test_prompt_hashes_match_protocol_7_sources(self) -> None:
        for name, expected in EXPECTED_PROMPT_HASHES.items():
            self.assertEqual(sha256_file(PACKAGE_ROOT / "prompts" / name), expected)

    def test_end_to_end_demo_matches_public_example(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            generated = Path(temporary) / "example"
            result = run_demo(generated)
            self.assertEqual(result["summary"]["frozen_claim_count"], 1)
            self.assertEqual(result["summary"]["frozen_evidence_count"], 1)
            public_example = PACKAGE_ROOT / "examples"
            expected_files = sorted(
                path.relative_to(generated)
                for path in generated.rglob("*")
                if path.is_file()
            )
            self.assertTrue(expected_files)
            for relative in expected_files:
                self.assertEqual(
                    (generated / relative).read_bytes(),
                    (public_example / relative).read_bytes(),
                    str(relative),
                )

    def test_modified_excerpt_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary) / "example"
            result = run_demo(base)
            bad_response = copy.deepcopy(result["classification_response"])
            bad_response["classifications"][0]["verified_excerpt"] = "text not in source"
            bad_path = base / "responses" / "bad_classification.jsonl"
            write_jsonl(bad_path, [bad_response])
            with self.assertRaises(AuditError):
                validate_responses(
                    result["classification_request_path"],
                    bad_path,
                    base / "expected" / "bad.jsonl",
                )

    def test_existing_frozen_interface(self) -> None:
        result = verify_reference_frozen_inputs(PACKAGE_ROOT.parent / "input")
        self.assertEqual(result["frozen_claim_count"], 498)
        self.assertEqual(result["frozen_evidence_count"], 2141)
        self.assertEqual(result["accepted_auditor_packets"], 114)
        self.assertAlmostEqual(
            result["semantic_discount_for_unreviewed_sources"],
            0.9674014658411585,
            places=15,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh-examples",
        action="store_true",
        help="Regenerate the deterministic, synthetic examples shipped with the package.",
    )
    args, remaining = parser.parse_known_args()
    if args.refresh_examples:
        run_demo(PACKAGE_ROOT / "examples")
        print(f"Refreshed examples in {PACKAGE_ROOT / 'examples'}")
        return
    unittest.main(argv=[__file__, *remaining])


if __name__ == "__main__":
    main()
