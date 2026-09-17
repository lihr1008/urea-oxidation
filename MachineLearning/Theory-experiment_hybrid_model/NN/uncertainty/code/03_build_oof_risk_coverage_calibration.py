from uncertainty_common import (
    build_oof_risk_coverage_calibration,
    load_config,
    write_grid_uncertainty_overview,
)


def main():
    _, config = load_config()
    result_dir = build_oof_risk_coverage_calibration(config)
    overview_path = write_grid_uncertainty_overview(config)
    print(f"Saved OOF risk-coverage calibration outputs to {result_dir}")
    print(f"Updated grid uncertainty overview at {overview_path}")


if __name__ == "__main__":
    main()
