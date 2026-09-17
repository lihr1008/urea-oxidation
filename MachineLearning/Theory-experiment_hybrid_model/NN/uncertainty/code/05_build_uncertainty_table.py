from uncertainty_common import build_uncertainty_table, load_config


def main():
    _, config = load_config()
    result_dir = build_uncertainty_table(config)
    print(f"Saved uncertainty table to {result_dir / 'uncertainty_table.xlsx'}")


if __name__ == "__main__":
    main()
