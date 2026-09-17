from uncertainty_common import (
    load_config,
    train_cv_seed_ensemble,
    write_grid_uncertainty_overview,
)


def main():
    _, config = load_config()
    result_dir = train_cv_seed_ensemble(config)
    overview_path = write_grid_uncertainty_overview(config)
    print(f"Saved repeated 10-fold CV seed ensemble outputs to {result_dir}")
    print(f"Updated grid uncertainty overview at {overview_path}")


if __name__ == "__main__":
    main()
