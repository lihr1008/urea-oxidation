from pathlib import Path

import pandas as pd


CODE_DIR = Path(__file__).resolve().parent
BASE_DIR = CODE_DIR.parent
RATIO_COLUMNS = ["Mn", "Fe", "Co", "Ni", "Zn"]
ZONE_ORDER = ["Z1", "Z2", "Z3", "Z4"]
ZONE_COLORS = {
    "Z1": "#54A24B",
    "Z2": "#B79A20",
    "Z3": "#E45756",
    "Z4": "#8E6C8A",
}
TOP_RANK_LIMIT = 3000
TOP_RANK_BIN_SIZE = 30

RESULT_DIR = BASE_DIR / "result"
OVERVIEW_PATH = RESULT_DIR / "grid_uncertainty_overview.xlsx"
OUTPUT_DIR = RESULT_DIR / "grid_uncertainty_count"


def build_rank_table(grid):
    required = RATIO_COLUMNS + ["predicted_potential", "Zone"]
    missing = [column for column in required if column not in grid.columns]
    if missing:
        raise ValueError(f"grid_uncertainty_overview.xlsx is missing columns: {missing}")

    ranking = grid[required].copy()
    ranking["_original_order"] = range(len(ranking))

    global_sorted = ranking.sort_values(
        ["predicted_potential", "_original_order"],
        ascending=[True, True],
        kind="mergesort",
    )
    ranking.loc[global_sorted.index, "global_rank"] = range(1, len(global_sorted) + 1)

    ranking["zone_rank"] = 0
    for zone in ZONE_ORDER:
        zone_sorted = ranking.loc[ranking["Zone"] == zone].sort_values(
            ["predicted_potential", "_original_order"],
            ascending=[True, True],
            kind="mergesort",
        )
        ranking.loc[zone_sorted.index, "zone_rank"] = range(1, len(zone_sorted) + 1)

    ranking["global_rank"] = ranking["global_rank"].astype(int)
    ranking["zone_rank"] = ranking["zone_rank"].astype(int)
    ranking = ranking.sort_values("global_rank", kind="mergesort")
    return ranking[RATIO_COLUMNS + ["predicted_potential", "global_rank", "zone_rank", "Zone"]]


def save_rank_workbook(ranking, output_path):
    with pd.ExcelWriter(output_path) as writer:
        ranking.to_excel(writer, sheet_name="all", index=False)
        for zone in ZONE_ORDER:
            zone_rank = ranking.loc[ranking["Zone"] == zone].sort_values(
                "zone_rank",
                kind="mergesort",
            )
            zone_rank.to_excel(writer, sheet_name=zone, index=False)


def build_origin_bin_table(ranking, rank_limit, bin_size):
    top = ranking.loc[ranking["global_rank"] <= int(rank_limit)].copy()
    if len(top) != int(rank_limit):
        raise ValueError(f"Expected {rank_limit} Top-ranked points, found {len(top)}.")
    if int(rank_limit) % int(bin_size) != 0:
        raise ValueError("TOP_RANK_LIMIT must be divisible by TOP_RANK_BIN_SIZE.")

    top["rank_bin"] = ((top["global_rank"] - 1) // int(bin_size) + 1).astype(int)
    counts = top.groupby(["rank_bin", "Zone"]).size().unstack(fill_value=0)
    bin_count = int(rank_limit) // int(bin_size)
    counts = counts.reindex(range(1, bin_count + 1), fill_value=0)
    counts = counts.reindex(columns=ZONE_ORDER, fill_value=0)

    output = pd.DataFrame({"rank_bin": range(1, bin_count + 1)})
    output["rank_start"] = (output["rank_bin"] - 1) * int(bin_size) + 1
    output["rank_end"] = output["rank_bin"] * int(bin_size)
    output["rank_center"] = (output["rank_start"] + output["rank_end"]) / 2.0
    output["rank_label"] = output["rank_start"].astype(str) + "-" + output["rank_end"].astype(str)
    output["total_count"] = counts.sum(axis=1).to_numpy(dtype=int)
    for zone in ZONE_ORDER:
        output[f"{zone}_count"] = counts[zone].to_numpy(dtype=int)
    for zone in ZONE_ORDER:
        output[f"{zone}_fraction"] = output[f"{zone}_count"] / output["total_count"]

    count_columns = [f"{zone}_count" for zone in ZONE_ORDER]
    fraction_columns = [f"{zone}_fraction" for zone in ZONE_ORDER]
    if int(output["total_count"].sum()) != int(rank_limit):
        raise ValueError("Origin bin table does not contain the requested number of ranked points.")
    if not (output["total_count"] == int(bin_size)).all():
        raise ValueError("Every Top-rank bin must contain exactly TOP_RANK_BIN_SIZE points.")
    if not (output[count_columns].sum(axis=1) == output["total_count"]).all():
        raise ValueError("Z1-Z4 counts do not sum to total_count in every bin.")
    if ((output[fraction_columns].sum(axis=1) - 1.0).abs() > 1e-12).any():
        raise ValueError("Z1-Z4 fractions do not sum to 1 in every bin.")
    return top, output


def save_origin_csv(bin_table, output_path):
    bin_table.to_csv(output_path, index=False, encoding="utf-8-sig", float_format="%.10g")


def save_rank_zone_distribution(plt, top_ranking, bin_table, output_path):
    fraction_columns = [f"{zone}_fraction" for zone in ZONE_ORDER]
    zone_counts = top_ranking["Zone"].value_counts().reindex(ZONE_ORDER, fill_value=0)
    labels = [f"{zone} (n={int(zone_counts[zone])})" for zone in ZONE_ORDER]

    fig, ax = plt.subplots(figsize=(13.5, 4.2), dpi=180)
    ax.stackplot(
        bin_table["rank_center"],
        *[bin_table[column].to_numpy(dtype=float) for column in fraction_columns],
        colors=[ZONE_COLORS[zone] for zone in ZONE_ORDER],
        labels=labels,
        step="mid",
        linewidth=0,
        rasterized=True,
    )
    ax.set_title(
        f"Zone distribution along global predicted-potential rank "
        f"(top-{TOP_RANK_LIMIT}, bin size = {TOP_RANK_BIN_SIZE})"
    )
    ax.set_xlabel("Global rank of predicted potential (ascending; 1 = lowest)")
    ax.set_ylabel("Zone proportion in each rank bin")
    ax.set_xlim(0.5, TOP_RANK_LIMIT + 0.5)
    ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.legend(frameon=False, loc="upper right", ncol=4)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)


def cleanup_legacy_outputs(current_png, current_csv):
    known_legacy_names = {
        "grid_uncertainty_count_bars.png",
        "grid_predicted_potential_global_rank_zone_distribution.png",
        "grid_ilr_pca_zone_with_experiments.png",
        "grid_predicted_potential_rank_scatter.png",
    }
    for name in known_legacy_names:
        path = OUTPUT_DIR / name
        if path.exists():
            path.unlink()
    for path in OUTPUT_DIR.glob("grid_predicted_potential_top*_rank_zone_distribution.png"):
        if path.resolve() != current_png.resolve():
            path.unlink()
    for path in OUTPUT_DIR.glob("grid_predicted_potential_top*_zone_distribution_origin.csv"):
        if path.resolve() != current_csv.resolve():
            path.unlink()


def main():
    if not OVERVIEW_PATH.exists():
        raise FileNotFoundError(f"Missing grid overview file: {OVERVIEW_PATH}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    grid = pd.read_excel(OVERVIEW_PATH, sheet_name="grid")
    ranking = build_rank_table(grid)
    top_ranking, origin_bins = build_origin_bin_table(
        ranking,
        TOP_RANK_LIMIT,
        TOP_RANK_BIN_SIZE,
    )

    rank_workbook_path = OUTPUT_DIR / "grid_predicted_potential_rankings.xlsx"
    top_png_path = OUTPUT_DIR / f"grid_predicted_potential_top{TOP_RANK_LIMIT}_rank_zone_distribution.png"
    origin_csv_path = OUTPUT_DIR / f"grid_predicted_potential_top{TOP_RANK_LIMIT}_zone_distribution_origin.csv"

    save_rank_workbook(ranking, rank_workbook_path)
    save_origin_csv(origin_bins, origin_csv_path)

    # Loading PyTorch first avoids a Windows runtime-library conflict in the
    # project Conda environment when Matplotlib calls NumPy linear algebra.
    try:
        import torch  # noqa: F401
    except ModuleNotFoundError:
        pass
    import matplotlib.pyplot as plt

    save_rank_zone_distribution(plt, top_ranking, origin_bins, top_png_path)
    cleanup_legacy_outputs(top_png_path, origin_csv_path)

    zone_counts = top_ranking["Zone"].value_counts().reindex(ZONE_ORDER, fill_value=0)
    print(f"Saved ranking workbook to {rank_workbook_path}")
    print(f"Saved Origin plot data to {origin_csv_path}")
    print(f"Saved top-{TOP_RANK_LIMIT} rank distribution to {top_png_path}")
    print(f"Top-{TOP_RANK_LIMIT} zone counts: {zone_counts.to_dict()}")


if __name__ == "__main__":
    main()
