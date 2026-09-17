from uncertainty_common import load_config, recommend_high_potential_candidates


def main():
    _, config = load_config()
    result_dir = recommend_high_potential_candidates(config)
    print(f"Saved high-potential candidate recommendations to {result_dir}")


if __name__ == "__main__":
    main()
