"""Recompute the knowledge-fusion and sensitivity result tables."""
from evidence_fusion import run_pipeline
from fusion_common import RESULT_DIR, verify_inputs


if __name__ == "__main__":
    verify_inputs(write_result=True)
    outputs = run_pipeline(write_results=True)
    print(
        "Knowledge fusion complete: "
        f"{len(outputs['run_candidate_decisions'])} nominal run decisions, "
        f"{len(outputs['run_dimension_masses'])} dimension rows, "
        f"{len(outputs['candidate_summary'])} candidates."
    )
    print(f"Machine-readable results: {RESULT_DIR.name}/")
