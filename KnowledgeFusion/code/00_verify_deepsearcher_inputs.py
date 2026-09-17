"""Verify the immutable DeepSearcher outputs and frozen audit ledgers."""
from fusion_common import verify_inputs


if __name__ == "__main__":
    summary = verify_inputs(write_result=True)
    print(
        "Verified "
        f"{summary['assisted_output_count']} assisted outputs, "
        f"{summary['control_output_count']} controls, and "
        f"{summary['corpus_partition_count']} corpus partitions."
    )
