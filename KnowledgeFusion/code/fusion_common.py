"""Shared paths, deterministic I/O, and fixed-input verification."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence


class FusionError(RuntimeError):
    """Raised when a frozen input or fusion invariant is violated."""


CODE_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = CODE_DIR.parent
REPO_ROOT = PACKAGE_DIR.parent
CONFIG_PATH = PACKAGE_DIR / "config" / "fusion_config.json"
INPUT_DIR = PACKAGE_DIR / "input"
RESULT_DIR = PACKAGE_DIR / "result"
MANIFEST_PATH = INPUT_DIR / "deepsearcher_manifest.json"
VERIFICATION_PATH = RESULT_DIR / "input_verification.json"


def text_of(value: Any) -> str:
    return "" if value is None else str(value).strip()


def as_bool(value: Any) -> bool:
    return text_of(value).lower() in {"1", "true", "yes", "y"}


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise FusionError(f"Expected a JSON object: {path.name}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    if not fields:
        raise FusionError(f"Refusing to write an empty table: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_text_sha256(path: Path) -> str:
    text = path.read_text(encoding="utf-8-sig")
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def safe_repo_path(relative_path: str) -> Path:
    if not relative_path or Path(relative_path).is_absolute() or ".." in Path(relative_path).parts:
        raise FusionError(f"Unsafe repository-relative path: {relative_path!r}")
    path = (REPO_ROOT / Path(relative_path)).resolve()
    try:
        path.relative_to(REPO_ROOT.resolve())
    except ValueError as exc:
        raise FusionError(f"Path escapes the repository: {relative_path}") from exc
    return path


def _expected_file_set(rows: Iterable[dict[str, Any]]) -> set[str]:
    return {Path(text_of(row.get("relative_path"))).as_posix() for row in rows}


def _actual_text_file_set(relative_root: str) -> set[str]:
    root = safe_repo_path(relative_root)
    return {
        path.relative_to(REPO_ROOT).as_posix()
        for path in root.rglob("*.txt")
        if path.is_file()
    }


def verify_inputs(write_result: bool = True) -> dict[str, Any]:
    """Verify all fixed DeepSearcher and frozen-audit inputs exactly."""
    if not MANIFEST_PATH.exists():
        raise FusionError("Missing input/deepsearcher_manifest.json")
    manifest = load_json(MANIFEST_PATH)
    runs = manifest.get("runs", [])
    corpus = manifest.get("corpus", [])
    if not isinstance(runs, list) or not isinstance(corpus, list):
        raise FusionError("Manifest runs and corpus entries must be lists")

    required_groups = {
        "with_deepsearcher": "Deepsearcher/result/elements",
        "without_deepsearcher": "Deepsearcher/without_deepsearcher",
    }
    group_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    record_ids: set[str] = set()
    for row in runs:
        record_id = text_of(row.get("record_id"))
        group = text_of(row.get("source_group"))
        model = text_of(row.get("model"))
        if not record_id or record_id in record_ids:
            raise FusionError(f"Duplicate or blank record_id in manifest: {record_id!r}")
        if group not in required_groups:
            raise FusionError(f"Unexpected source_group for {record_id}: {group}")
        if model not in {"deepseek", "gemini", "gpt", "qwen"}:
            raise FusionError(f"Unexpected model for {record_id}: {model}")
        record_ids.add(record_id)
        group_rows[group].append(row)
        path = safe_repo_path(text_of(row.get("relative_path")))
        if not path.is_file():
            raise FusionError(f"Missing fixed output: {row.get('relative_path')}")
        if sha256_file(path) != text_of(row.get("sha256")):
            raise FusionError(f"SHA-256 mismatch: {row.get('relative_path')}")

    model_distribution: dict[str, dict[str, int]] = {}
    expected_counts = Counter({"deepseek": 10, "gemini": 10, "gpt": 10, "qwen": 10})
    for group, relative_root in required_groups.items():
        rows = group_rows[group]
        expected = _expected_file_set(rows)
        actual = _actual_text_file_set(relative_root)
        if expected != actual:
            missing = sorted(expected - actual)
            additional = sorted(actual - expected)
            raise FusionError(
                f"Fixed output set mismatch for {group}; missing={missing}, additional={additional}"
            )
        counts = Counter(text_of(row.get("model")) for row in rows)
        if len(rows) != 40 or counts != expected_counts:
            raise FusionError(f"Expected 40 files and 10 runs per model for {group}; found {dict(counts)}")
        model_distribution[group] = dict(sorted(counts.items()))

    expected_corpus = _expected_file_set(corpus)
    actual_corpus = _actual_text_file_set("Deepsearcher/data")
    if len(corpus) != 10 or expected_corpus != actual_corpus:
        raise FusionError("Expected exactly the 10 manifested Deepsearcher/data text partitions")
    for row in corpus:
        path = safe_repo_path(text_of(row.get("relative_path")))
        if normalized_text_sha256(path) != text_of(row.get("normalized_text_sha256")):
            raise FusionError(f"Normalized corpus hash mismatch: {row.get('relative_path')}")

    frozen_bindings = manifest.get("frozen_inputs", {})
    expected_frozen = {
        "frozen_claims.csv",
        "frozen_evidence.csv",
        "frozen_review_summary.json",
        "fusion_config.json",
    }
    if set(frozen_bindings) != expected_frozen:
        raise FusionError("Manifest must bind all frozen ledgers and the fusion configuration")
    frozen_paths = {
        "frozen_claims.csv": INPUT_DIR / "frozen_claims.csv",
        "frozen_evidence.csv": INPUT_DIR / "frozen_evidence.csv",
        "frozen_review_summary.json": INPUT_DIR / "frozen_review_summary.json",
        "fusion_config.json": CONFIG_PATH,
    }
    for name, path in frozen_paths.items():
        if not path.is_file() or sha256_file(path) != text_of(frozen_bindings.get(name)):
            raise FusionError(f"Frozen-input binding mismatch: {name}")

    result = {
        "status": "verified",
        "manifest_sha256": sha256_file(MANIFEST_PATH),
        "assisted_output_count": len(group_rows["with_deepsearcher"]),
        "control_output_count": len(group_rows["without_deepsearcher"]),
        "corpus_partition_count": len(corpus),
        "model_distribution": model_distribution,
        "corpus_hash_mode": "UTF-8 text with CRLF and CR normalized to LF",
        "frozen_input_count": len(frozen_bindings),
        "all_hashes_match": True,
        "no_fallback_used": True,
    }
    if write_result:
        write_json(VERIFICATION_PATH, result)
    return result


def load_config() -> dict[str, Any]:
    return load_json(CONFIG_PATH)
