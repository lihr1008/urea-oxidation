"""Deterministic regression checks for the lightweight fusion package."""
from __future__ import annotations

from evidence_fusion import run_pipeline
from fusion_common import FusionError, RESULT_DIR, load_config, read_csv, verify_inputs


def assert_close(left: object, right: object, tolerance: float = 1e-12) -> None:
    if left == right:
        return
    try:
        if abs(float(left) - float(right)) <= tolerance:
            return
    except (TypeError, ValueError):
        pass
    raise AssertionError(f"Values differ: {left!r} != {right!r}")


def compare_saved_table(name: str, generated: list[dict[str, object]]) -> None:
    saved = read_csv(RESULT_DIR / f"{name}.csv")
    if len(saved) != len(generated):
        raise AssertionError(f"{name}: row count changed")
    for index, (expected, actual) in enumerate(zip(generated, saved)):
        if set(expected) != set(actual):
            raise AssertionError(f"{name} row {index}: columns changed")
        for key in expected:
            assert_close(expected[key], actual[key])


def main() -> None:
    verified = verify_inputs(write_result=False)
    if not verified["all_hashes_match"]:
        raise AssertionError("Input verification did not pass")
    outputs = run_pipeline(write_results=False)
    config = load_config()
    expected = config["reference_expectations"]

    if len(outputs["run_candidate_decisions"]) != expected["assisted_run_count"]:
        raise AssertionError("Expected 40 nominal run decisions")
    if len(outputs["run_dimension_masses"]) != expected["nominal_dimension_count"]:
        raise AssertionError("Expected 160 nominal dimension masses")
    if len(outputs["parameter_sensitivity"]) != expected["scenario_count"] * expected["candidate_count"]:
        raise AssertionError("Expected five candidate summaries in each of 11 scenarios")

    for row in outputs["run_dimension_masses"]:
        total = float(row["R"]) + float(row["NR"]) + float(row["U"])
        if abs(total - 1.0) > 1e-12:
            raise AssertionError("R + NR + U differs from one")

    priority = next(
        row
        for row in outputs["candidate_summary"]
        if row["candidate"] == expected["priority_candidate"]
    )
    if priority["generation_count"] != expected["priority_generation_count"]:
        raise AssertionError("Priority candidate generation count changed")
    if priority["supported_run_count"] != expected["priority_supported_run_count"]:
        raise AssertionError("Priority candidate gate-pass count changed")

    for name, rows in outputs.items():
        compare_saved_table(name, rows)
    print("All knowledge-fusion verification checks passed.")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, FusionError) as exc:
        raise SystemExit(str(exc)) from exc
