from pathlib import Path
import sys

from nn_common import (
    create_analysis_from_space_prediction,
    get_run_dirs_by_id,
    run_ilr_pca,
)

RUN_ID = "experiment_only_reference"


def main(run_id=RUN_ID):
    base_dir = Path(__file__).resolve().parent
    _, _, result_dir = get_run_dirs_by_id(base_dir, run_id)
    analysis_path = result_dir / "analysis.xlsx"
    create_analysis_from_space_prediction(
        result_dir / "space_one_percent.xlsx", analysis_path
    )
    run_ilr_pca(analysis_path, result_dir)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else RUN_ID)
