from __future__ import annotations

import numpy as np
import pandas as pd

from run_shap_analysis import (
    DESCRIPTOR_DISPLAY_COLUMNS,
    DESCRIPTOR_DISPLAY_TO_PLOT,
    DESCRIPTOR_PHYSICAL_COLUMNS,
    DESCRIPTOR_PLOT_LABELS,
    RATIO_COLUMNS,
    RESULT_DIR,
    importance_frame,
    save_combined_importance_plot,
    save_summary_plot,
    set_plot_style,
)


SOURCE_WORKBOOK = RESULT_DIR / "descriptor_to_E10_shap_data.xlsx"
IMPORTANCE_OUTPUT = (
    RESULT_DIR / "descriptor_to_E10_shap_all_excluding_monometallic_Mn.png"
)
SUMMARY_OUTPUT = (
    RESULT_DIR
    / "descriptor_to_E10_shap_summary_plot_excluding_monometallic_Mn.png"
)
EXPECTED_REMOVED_SOURCE_INDEX = 105
EXPECTED_RETAINED_SAMPLES = 21


def require_columns(frame: pd.DataFrame, columns: list[str], sheet_name: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"Sheet {sheet_name!r} is missing columns: {missing}")


def monometallic_mn_mask(metals: pd.DataFrame) -> np.ndarray:
    other_metals = [metal for metal in RATIO_COLUMNS if metal != "Mn"]
    return np.isclose(metals["Mn"].to_numpy(dtype=float), 1.0) & np.isclose(
        metals[other_metals].to_numpy(dtype=float), 0.0
    ).all(axis=1)


def main() -> None:
    metals = pd.read_excel(SOURCE_WORKBOOK, sheet_name="metals")
    shap_frame = pd.read_excel(SOURCE_WORKBOOK, sheet_name="shap")
    descriptor_frame = pd.read_excel(
        SOURCE_WORKBOOK,
        sheet_name="x_val_physical_raw",
    )

    require_columns(metals, ["source_index", *RATIO_COLUMNS], "metals")
    require_columns(shap_frame, DESCRIPTOR_DISPLAY_COLUMNS, "shap")
    require_columns(
        descriptor_frame,
        DESCRIPTOR_PHYSICAL_COLUMNS,
        "x_val_physical_raw",
    )

    row_counts = {len(metals), len(shap_frame), len(descriptor_frame)}
    if len(row_counts) != 1:
        raise ValueError(
            "The metals, SHAP, and descriptor sheets do not have aligned row counts."
        )

    remove_mask = monometallic_mn_mask(metals)
    removed_rows = metals.loc[remove_mask]
    if len(removed_rows) != 1:
        raise ValueError(
            "Expected exactly one monometallic Mn sample, "
            f"but found {len(removed_rows)}."
        )

    removed_source_index = int(removed_rows.iloc[0]["source_index"])
    if removed_source_index != EXPECTED_REMOVED_SOURCE_INDEX:
        raise ValueError(
            "The monometallic Mn sample has unexpected source_index "
            f"{removed_source_index}; expected {EXPECTED_REMOVED_SOURCE_INDEX}."
        )

    keep_mask = ~remove_mask
    retained_count = int(keep_mask.sum())
    if retained_count != EXPECTED_RETAINED_SAMPLES:
        raise ValueError(
            f"Expected {EXPECTED_RETAINED_SAMPLES} retained samples, "
            f"but found {retained_count}."
        )

    shap_values = shap_frame[DESCRIPTOR_DISPLAY_COLUMNS].to_numpy(dtype=float)[
        keep_mask
    ]
    descriptor_values = descriptor_frame[
        DESCRIPTOR_PHYSICAL_COLUMNS
    ].to_numpy(dtype=float)[keep_mask]

    importance = importance_frame(shap_values, DESCRIPTOR_DISPLAY_COLUMNS)
    plot_importance = importance.copy()
    plot_importance["feature"] = plot_importance["feature"].map(
        DESCRIPTOR_DISPLAY_TO_PLOT
    )

    set_plot_style()
    save_summary_plot(
        shap_values,
        descriptor_values,
        DESCRIPTOR_PLOT_LABELS,
        SUMMARY_OUTPUT,
    )
    save_combined_importance_plot(
        plot_importance,
        "Mean |SHAP Value|",
        IMPORTANCE_OUTPUT,
    )

    print(f"Removed monometallic Mn source_index: {removed_source_index}")
    print(f"Retained samples: {retained_count}")
    print(importance.to_string(index=False))
    print(f"Saved: {IMPORTANCE_OUTPUT}")
    print(f"Saved: {SUMMARY_OUTPUT}")


if __name__ == "__main__":
    main()
