from uncertainty_common import (
    build_applicability_domain,
    load_config,
    write_grid_uncertainty_overview,
)


def main():
    _, config = load_config()
    result_dir = build_applicability_domain(config)
    overview_path = write_grid_uncertainty_overview(config)
    print(f"Saved applicability-domain outputs to {result_dir}")
    print(f"Updated grid uncertainty overview at {overview_path}")


if __name__ == "__main__":
    main()
