from uncertainty_common import (
    build_descriptor_error_propagation,
    load_config,
)


def main():
    _, config = load_config()
    result_dir = build_descriptor_error_propagation(config)
    print(f"Saved descriptor propagation analysis to {result_dir / 'descriptor_error_propagation.xlsx'}")


if __name__ == "__main__":
    main()
