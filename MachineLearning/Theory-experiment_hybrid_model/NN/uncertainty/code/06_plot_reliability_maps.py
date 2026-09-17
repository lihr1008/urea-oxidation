from uncertainty_common import load_config, plot_reliability_2x4


def main():
    _, config = load_config()
    result_dir = plot_reliability_2x4(config)
    print(f"Saved reliability map to {result_dir / 'reliability_2x4_map.png'}")


if __name__ == "__main__":
    main()
