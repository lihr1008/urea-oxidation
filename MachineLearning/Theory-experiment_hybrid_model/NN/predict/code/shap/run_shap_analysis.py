from __future__ import annotations

import logging
import os
import pickle
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from matplotlib.cm import ScalarMappable
from sklearn.metrics import pairwise_distances_argmin_min
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from torch.nn import functional as F


def add_conda_dll_directories():
    if os.name != "nt":
        return
    env_dir = Path(sys.executable).resolve().parent
    for path in [env_dir, env_dir / "Library" / "bin", env_dir / "Scripts"]:
        if not path.exists():
            continue
        os.environ["PATH"] = str(path) + os.pathsep + os.environ.get("PATH", "")
        add_dll_directory = getattr(os, "add_dll_directory", None)
        if add_dll_directory is not None:
            try:
                add_dll_directory(str(path))
            except OSError:
                pass


add_conda_dll_directories()

import shap


logging.getLogger("shap").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message="Trying to unpickle estimator StandardScaler.*")
warnings.filterwarnings("ignore", message="Using .* background data samples could cause slower run times.*")


SHAP_DIR = Path(__file__).resolve().parent
CODE_DIR = SHAP_DIR.parent
PREDICT_DIR = CODE_DIR.parent
NN_DIR = PREDICT_DIR.parent
ML_ROOT = NN_DIR.parent.parent
REPO_ROOT = ML_ROOT.parent
RESULT_DIR = SHAP_DIR / "results"

WEIGHT_DIR = PREDICT_DIR / "results"
EXPERIMENT_XLSX = ML_ROOT / "data" / "experiment.xlsx"
THEORY_XLSX = ML_ROOT / "data" / "metals_activity.xlsx"

THEORY_SPLIT_SEED = 0
EXPERIMENT_SPLIT_SEED = 106
TEST_SIZE = 0.2
GLOBAL_BACKGROUND_SIZES = (100, 200, 500)
GLOBAL_REFERENCE_BACKGROUND_SIZE = 500
SHAP_CHUNK_SIZE = 100
GLOBAL_RANK_CORRELATION_THRESHOLD = 0.90
GLOBAL_IMPORTANCE_PERCENT_TOLERANCE = 2.0
GLOBAL_BASELINE_STD_TOLERANCE = 0.01
GLOBAL_MAJOR_IMPORTANCE_PERCENT = 5.0

RATIO_COLUMNS = ["Mn", "Fe", "Co", "Ni", "Zn"]
DESCRIPTOR_MODEL_COLUMNS = ["ade_urea", "ade_CO", "e_urea", "e_CO"]
DESCRIPTOR_SCALER_COLUMNS = ["\u0394E_urea", "\u0394E_CO", "\u0394e_urea", "\u0394q_CO"]
DESCRIPTOR_PHYSICAL_COLUMNS = ["\u0394E_urea", "\u0394E_CO", "\u0394q_urea", "\u0394q_CO"]
DESCRIPTOR_DISPLAY_COLUMNS = [
    "ade_urea / \u0394E_urea",
    "ade_CO / \u0394E_CO",
    "e_urea / \u0394q_urea",
    "e_CO / \u0394q_CO",
]
DESCRIPTOR_PLOT_LABELS = [
    r"$\Delta E_{\mathrm{*urea}}$",
    r"$\Delta E_{\mathrm{*CO}}$",
    r"$\Delta q_{\mathrm{*urea}}$",
    r"$\Delta q_{\mathrm{*CO}}$",
]
DESCRIPTOR_DISPLAY_TO_PLOT = dict(zip(DESCRIPTOR_DISPLAY_COLUMNS, DESCRIPTOR_PLOT_LABELS))
TARGETS = [
    {"model_column": "ade_urea", "physical_column": "\u0394E_urea", "slug": "ade_urea", "index": 0},
    {"model_column": "ade_CO", "physical_column": "\u0394E_CO", "slug": "ade_CO", "index": 1},
    {"model_column": "e_urea", "physical_column": "\u0394q_urea", "slug": "e_urea", "index": 2},
    {"model_column": "e_CO", "physical_column": "\u0394q_CO", "slug": "e_CO", "index": 3},
]


class PreModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.Layer1 = nn.Linear(in_features=5, out_features=512)
        self.Layer2 = nn.Linear(in_features=512, out_features=512)
        self.Layer3 = nn.Linear(in_features=512, out_features=4)

    def forward(self, x):
        x = F.relu(self.Layer1(x))
        x = F.relu(self.Layer2(x))
        return self.Layer3(x)


class DescriptorToExperimentModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.Layer1 = nn.Linear(in_features=4, out_features=512)
        self.Layer2 = nn.Linear(in_features=512, out_features=512)
        self.Layer3 = nn.Linear(in_features=512, out_features=512)
        self.Layer4 = nn.Linear(in_features=512, out_features=1)
        self.dropout = nn.Dropout(0.4)

    def forward(self, x):
        x = self.dropout(F.relu(self.Layer1(x)))
        x = self.dropout(F.relu(self.Layer2(x)))
        x = self.dropout(F.relu(self.Layer3(x)))
        return self.Layer4(x)


def load_pickle(path: Path):
    with open(path, "rb") as f:
        return pickle.load(f)


def repository_relative(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def load_assets():
    premodel = PreModel()
    premodel.load_state_dict(torch.load(WEIGHT_DIR / "PBA-ml-pretrain_model_weights.pth", map_location="cpu"))
    premodel.eval()

    model = DescriptorToExperimentModel()
    model.load_state_dict(torch.load(WEIGHT_DIR / "train_model_weights.pth", map_location="cpu"))
    model.eval()

    for network in [premodel, model]:
        for param in network.parameters():
            param.requires_grad = False

    return {
        "premodel": premodel,
        "model": model,
        "metal_scaler": load_pickle(WEIGHT_DIR / "pretrain_norm_x.pkl"),
        "descriptor_scaler": load_pickle(WEIGHT_DIR / "pretrain_norm_y.pkl"),
        "target_scaler": load_pickle(WEIGHT_DIR / "train_norm_y.pkl"),
    }


def read_metal_frame(path: Path, sheet_name: str) -> pd.DataFrame:
    metals = pd.read_excel(path, sheet_name=sheet_name)
    missing = [col for col in RATIO_COLUMNS if col not in metals.columns]
    if missing:
        raise ValueError(f"{path} {sheet_name} sheet is missing columns: {missing}")
    metals = metals[RATIO_COLUMNS].copy()
    metals.insert(0, "source_index", metals.index.to_numpy(dtype=int))
    return metals


def split_metal_frame(metals: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    train, test = train_test_split(metals, test_size=TEST_SIZE, random_state=seed)
    return train.reset_index(drop=True), test.reset_index(drop=True)


def load_theory_split(assets: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    metals = read_metal_frame(THEORY_XLSX, "metals")
    train, test = split_metal_frame(metals, THEORY_SPLIT_SEED)
    fitted_scaler = StandardScaler().fit(train[RATIO_COLUMNS])
    loaded_scaler = assets["metal_scaler"]
    mean_error = float(np.max(np.abs(fitted_scaler.mean_ - loaded_scaler.mean_)))
    scale_error = float(np.max(np.abs(fitted_scaler.scale_ - loaded_scaler.scale_)))
    if not np.allclose(fitted_scaler.mean_, loaded_scaler.mean_, rtol=0.0, atol=1e-12):
        raise ValueError("The theory split mean does not match the fixed metal scaler.")
    if not np.allclose(fitted_scaler.scale_, loaded_scaler.scale_, rtol=0.0, atol=1e-12):
        raise ValueError("The theory split scale does not match the fixed metal scaler.")
    metadata = {
        "theory_rows_total": int(len(metals)),
        "theory_train_count": int(len(train)),
        "theory_test_count": int(len(test)),
        "theory_split_seed": THEORY_SPLIT_SEED,
        "theory_test_size": TEST_SIZE,
        "metal_scaler_mean_max_abs_error": mean_error,
        "metal_scaler_scale_max_abs_error": scale_error,
    }
    return train, test, metadata


def load_experiment_split() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    metals = read_metal_frame(EXPERIMENT_XLSX, "metals")
    train, test = split_metal_frame(metals, EXPERIMENT_SPLIT_SEED)
    metadata = {
        "experiment_rows_total": int(len(metals)),
        "experiment_train_count": int(len(train)),
        "experiment_test_count": int(len(test)),
        "experiment_split_seed": EXPERIMENT_SPLIT_SEED,
        "experiment_test_size": TEST_SIZE,
        "single_element_rows_total": int(((metals[RATIO_COLUMNS] > 0).sum(axis=1) <= 1).sum()),
        "single_element_rows_excluded": 0,
    }
    return train, test, metadata


def predict_descriptors_physical(metals, assets: dict) -> np.ndarray:
    metal_frame = pd.DataFrame(metals, columns=RATIO_COLUMNS)
    metal_scaled = assets["metal_scaler"].transform(metal_frame)
    with torch.no_grad():
        descriptor_scaled = assets["premodel"](torch.tensor(metal_scaled).to(torch.float32)).detach().numpy()
    return assets["descriptor_scaler"].inverse_transform(descriptor_scaled)


def predict_e10_from_descriptors(descriptors, assets: dict) -> np.ndarray:
    descriptor_frame = pd.DataFrame(descriptors, columns=DESCRIPTOR_SCALER_COLUMNS)
    descriptor_scaled = assets["descriptor_scaler"].transform(descriptor_frame)
    with torch.no_grad():
        pred_scaled = assets["model"](torch.tensor(descriptor_scaled).to(torch.float32)).detach().numpy()
    return assets["target_scaler"].inverse_transform(pred_scaled)


def normalize_shap_values(values) -> np.ndarray:
    if isinstance(values, list):
        values = values[0]
    values = np.asarray(values)
    if values.ndim == 3 and values.shape[-1] == 1:
        values = np.squeeze(values, axis=-1)
    if values.ndim != 2:
        raise ValueError(f"Expected 2-D SHAP values, got shape {values.shape}")
    return values


def normalize_multioutput_shap_values(values, n_outputs: int) -> np.ndarray:
    if isinstance(values, list):
        values = np.stack([np.asarray(value) for value in values], axis=-1)
    values = np.asarray(values)
    if values.ndim == 2 and n_outputs == 1:
        values = values[:, :, np.newaxis]
    if values.ndim != 3 or values.shape[-1] != n_outputs:
        raise ValueError(
            f"Expected SHAP values with shape (samples, features, {n_outputs}), "
            f"got {values.shape}"
        )
    return values


def exact_kernel_nsamples(n_features: int) -> int:
    return 2**n_features - 2


def compute_kernel_shap(
    model_fn,
    background,
    x_explain: np.ndarray,
    n_outputs: int,
    chunk_size: int = SHAP_CHUNK_SIZE,
) -> tuple[np.ndarray, np.ndarray]:
    explainer = shap.KernelExplainer(model_fn, background)
    chunks = []
    n_features = int(x_explain.shape[1])
    for start in range(0, len(x_explain), chunk_size):
        stop = min(start + chunk_size, len(x_explain))
        chunk_values = explainer.shap_values(
            x_explain[start:stop],
            nsamples=exact_kernel_nsamples(n_features),
            l1_reg=f"num_features({n_features})",
            silent=True,
            gc_collect=False,
        )
        chunks.append(normalize_multioutput_shap_values(chunk_values, n_outputs))
    expected_value = np.atleast_1d(np.asarray(explainer.expected_value, dtype=float))
    if expected_value.size != n_outputs:
        raise ValueError(
            f"Expected {n_outputs} SHAP base values, got {expected_value.shape}"
        )
    return np.concatenate(chunks, axis=0), expected_value


def validate_additivity(
    model_fn,
    x_explain: np.ndarray,
    shap_values: np.ndarray,
    expected_value: np.ndarray,
    output_names: list[str],
) -> pd.DataFrame:
    predictions = np.asarray(model_fn(x_explain), dtype=float)
    if predictions.ndim == 1:
        predictions = predictions[:, np.newaxis]
    reconstructed = expected_value.reshape(1, -1) + shap_values.sum(axis=1)
    errors = np.abs(predictions - reconstructed)
    rows = []
    for output_index, output_name in enumerate(output_names):
        passed = bool(
            np.allclose(
                predictions[:, output_index],
                reconstructed[:, output_index],
                rtol=1e-4,
                atol=1e-6,
            )
        )
        rows.append(
            {
                "output": output_name,
                "expected_value": float(expected_value[output_index]),
                "max_abs_error": float(errors[:, output_index].max()),
                "mean_abs_error": float(errors[:, output_index].mean()),
                "rtol": 1e-4,
                "atol": 1e-6,
                "passed": passed,
            }
        )
        if not passed:
            raise ValueError(f"SHAP additivity check failed for {output_name}.")
    return pd.DataFrame(rows)


def set_plot_style():
    plt.rcParams["font.serif"] = ["Arial"]
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams["axes.titleweight"] = "bold"
    plt.rcParams["xtick.labelsize"] = 12
    plt.rcParams["ytick.labelsize"] = 12
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.minor.width"] = 2.0
    plt.rcParams["ytick.minor.width"] = 1.0
    plt.rcParams["axes.unicode_minus"] = False


def importance_frame(shap_values: np.ndarray, features: list[str]) -> pd.DataFrame:
    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    total = float(mean_abs_shap.sum())
    fractions = mean_abs_shap / total if total else np.zeros_like(mean_abs_shap)
    return pd.DataFrame(
        {
            "feature": features,
            "mean_abs_shap": mean_abs_shap,
            "fraction": fractions,
            "percent": fractions * 100,
        }
    ).sort_values("mean_abs_shap", ascending=False, ignore_index=True)


def save_summary_plot(shap_values: np.ndarray, x_values: np.ndarray, feature_names: list[str], path: Path):
    plt.figure()
    shap.summary_plot(shap_values, x_values, feature_names=feature_names, show=False)
    plt.gca().tick_params(axis="y", pad=2)
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def save_combined_importance_plot(
    importance: pd.DataFrame,
    xlabel: str,
    path: Path,
    signed_effect: pd.Series | None = None,
):
    feature_names = importance["feature"].tolist()
    shap_series = importance["mean_abs_shap"].to_numpy(dtype=float)
    sorted_weights = importance["fraction"].to_numpy(dtype=float)
    num_vars = len(feature_names)
    signed_values = None
    if signed_effect is not None:
        signed_values = np.asarray([signed_effect.loc[feature] for feature in feature_names], dtype=float)

    pink = (253 / 255, 0 / 255, 89 / 255)
    blue = (0 / 255, 134 / 255, 250 / 255)
    cmap = mcolors.LinearSegmentedColormap.from_list("blue_pink_manual", [blue, pink], N=256)
    color_norm = plt.Normalize(shap_series.min(), shap_series.max())
    colors = cmap(color_norm(shap_series))

    fig = plt.figure(figsize=(16, 8))
    left_margin, right_margin, bottom_margin, top_margin = 0.08, 0.08, 0.12, 0.12
    space_between_plots, colorbar_width = 0.04, 0.02
    plot_bottom, plot_height = bottom_margin, 1.0 - bottom_margin - top_margin
    cbar_left = left_margin
    main_ax_left = cbar_left + colorbar_width + space_between_plots
    main_ax_width = 1.0 - main_ax_left - right_margin

    ax_cbar = fig.add_axes([cbar_left, plot_bottom, colorbar_width, plot_height])
    ax_bar = fig.add_axes([main_ax_left, plot_bottom, main_ax_width, plot_height])

    sm = ScalarMappable(cmap=cmap, norm=color_norm)
    cbar = fig.colorbar(sm, cax=ax_cbar, orientation="vertical")
    cbar.set_ticks([])
    cbar.ax.yaxis.set_ticks_position("left")
    ax_cbar.text(0.5, 1.01, "High", transform=ax_cbar.transAxes, ha="center", va="bottom", fontsize=24, fontweight="bold")
    ax_cbar.text(0.5, -0.01, "Low", transform=ax_cbar.transAxes, ha="center", va="top", fontsize=24, fontweight="bold")
    cbar.outline.set_visible(False)
    ax_cbar.text(-1.2, 0.5, "Mean |SHAP Value|", transform=ax_cbar.transAxes, fontsize=24, rotation=90, va="center", fontweight="bold")

    ax_bar.xaxis.tick_bottom()
    ax_bar.xaxis.set_label_position("bottom")
    ax_bar.invert_xaxis()
    ax_bar.barh(range(num_vars), shap_series, color=colors, height=0.6)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel(xlabel, size=24, labelpad=8, fontweight="bold")
    ax_bar.set_yticks([])
    ax_bar.spines[["left", "top"]].set_visible(False)
    ax_bar.spines["right"].set_position(("data", 0))
    ax_bar.spines["right"].set_visible(True)
    ax_bar.spines["bottom"].set_visible(True)
    ax_bar.tick_params(axis="x", which="major", direction="out", labelsize=20, length=6, pad=8)
    ax_bar.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax_bar.tick_params(axis="x", which="minor", direction="out", length=4)
    for label in ax_bar.get_xticklabels():
        label.set_fontweight("bold")

    x_text = -0.005 * shap_series.max()
    value_text_x = shap_series.max() * 0.04
    for i, feature in enumerate(feature_names[:num_vars]):
        if signed_values is None:
            label_text = feature
        else:
            direction = "\u2191" if signed_values[i] > 0 else "\u2193" if signed_values[i] < 0 else "\u2192"
            label_text = f"{feature} {direction}"
        ax_bar.text(x_text, i, label_text, ha="left", va="center", color="black", fontsize=24, fontweight="bold")
        if signed_values is not None:
            ax_bar.text(
                value_text_x,
                i,
                f"{signed_values[i]:+.4g}",
                ha="right",
                va="center",
                color="black",
                fontsize=18,
                fontweight="bold",
            )

    inset_left = main_ax_left - 0.2
    inset_bottom = plot_bottom - 0.08
    inset_size = min(main_ax_width, plot_height) * 0.75
    ax_radial_inset = fig.add_axes([inset_left, inset_bottom, inset_size, inset_size], projection="polar")
    ax_radial_inset.patch.set_alpha(0)

    percentages = sorted_weights * 100
    widths = sorted_weights * 2 * np.pi
    base_length, fixed_increment, colored_ring_width = 3.0, 0.5, 2.0
    total_lengths = [base_length + i * fixed_increment for i in range(num_vars)]
    inner_heights = [max(0, tl - colored_ring_width) for tl in total_lengths]
    inner_colors = (["#EAEAEA", "#FFFFFF"] * (num_vars // 2 + 1))[:num_vars]
    one_oclock_offset = np.pi / 21
    thetas = np.cumsum([0] + widths[:-1].tolist()) - one_oclock_offset

    ax_radial_inset.bar(thetas, inner_heights, width=widths, color=inner_colors, align="edge", edgecolor="white", linewidth=1.5)
    ax_radial_inset.bar(thetas, [colored_ring_width] * num_vars, width=widths, bottom=inner_heights, color=colors, align="edge", edgecolor="white", linewidth=1.5)
    for i in range(num_vars):
        label_angle_rad = thetas[i] + widths[i] / 1.8
        label_radius = total_lengths[i] + 1.2
        ax_radial_inset.text(label_angle_rad, label_radius, f"{percentages[i]:.1f}%", ha="center", va="center", fontsize=16, fontweight="bold")

    ax_radial_inset.set_yticklabels([])
    ax_radial_inset.set_xticklabels([])
    ax_radial_inset.spines["polar"].set_visible(False)
    ax_radial_inset.grid(False)
    ax_radial_inset.set_theta_zero_location("N")
    ax_radial_inset.set_theta_direction(-1)
    ax_radial_inset.set_ylim(0, max(total_lengths) + 3)

    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_matrix_heatmap(
    matrix: pd.DataFrame,
    path: Path,
    colorbar_label: str,
    signed: bool = False,
):
    values = matrix.loc[RATIO_COLUMNS, DESCRIPTOR_PHYSICAL_COLUMNS].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(10, 6))
    if signed:
        limit = float(np.max(np.abs(values)))
        limit = limit if limit > 0 else 1.0
        image = ax.imshow(values, cmap="coolwarm", aspect="auto", vmin=-limit, vmax=limit)
    else:
        image = ax.imshow(values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(DESCRIPTOR_PLOT_LABELS)), labels=DESCRIPTOR_PLOT_LABELS)
    ax.set_yticks(range(len(RATIO_COLUMNS)), labels=RATIO_COLUMNS)
    ax.tick_params(axis="x", pad=8)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            normalized = image.norm(values[row, column])
            text_color = "white" if normalized < 0.25 or normalized > 0.75 else "black"
            ax.text(
                column,
                row,
                f"{values[row, column]:.3g}",
                ha="center",
                va="center",
                color=text_color,
                fontsize=11,
                fontweight="bold",
            )
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label(colorbar_label, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_descriptor_heatmaps(
    importance_matrix: pd.DataFrame,
    signed_effect_matrix: pd.DataFrame,
    prefix: str = "",
):
    save_matrix_heatmap(
        importance_matrix,
        RESULT_DIR / f"{prefix}metal_to_descriptor_mean_abs_shap_heatmap.png",
        "Mean |SHAP Value|",
    )
    save_matrix_heatmap(
        signed_effect_matrix,
        RESULT_DIR / f"{prefix}metal_to_descriptor_signed_shap_heatmap.png",
        "Signed mean |SHAP Value|",
        signed=True,
    )


def background_quality_frames(
    theory_train: pd.DataFrame,
    background,
    metal_scaler,
    background_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    centers = np.asarray(background.data, dtype=float)
    weights = np.asarray(background.weights, dtype=float)
    if not np.isclose(weights.sum(), 1.0, rtol=0.0, atol=1e-12):
        raise ValueError(f"Background weights for K={background_size} do not sum to one.")
    composition_error = np.abs(centers.sum(axis=1) - 1.0)
    if composition_error.max() > 1e-10 or centers.min() < -1e-12:
        raise ValueError(f"Background centers for K={background_size} are invalid compositions.")

    center_frame = pd.DataFrame(centers, columns=RATIO_COLUMNS)
    center_frame.insert(0, "cluster_id", np.arange(len(center_frame), dtype=int))
    center_frame["weight"] = weights

    train_values = theory_train[RATIO_COLUMNS].to_numpy(dtype=float)
    train_mean = train_values.mean(axis=0)
    train_std = train_values.std(axis=0)
    background_mean = np.average(centers, axis=0, weights=weights)
    background_variance = np.average(
        (centers - background_mean.reshape(1, -1)) ** 2,
        axis=0,
        weights=weights,
    )
    background_std = np.sqrt(background_variance)
    moment_rows = []
    for index, metal in enumerate(RATIO_COLUMNS):
        moment_rows.append(
            {
                "K": background_size,
                "metal": metal,
                "train_mean": float(train_mean[index]),
                "background_weighted_mean": float(background_mean[index]),
                "mean_abs_error": float(abs(background_mean[index] - train_mean[index])),
                "train_std": float(train_std[index]),
                "background_weighted_std": float(background_std[index]),
                "std_relative_error": float(
                    abs(background_std[index] - train_std[index])
                    / max(train_std[index], 1e-12)
                ),
            }
        )

    train_scaled = metal_scaler.transform(theory_train[RATIO_COLUMNS])
    center_scaled = metal_scaler.transform(pd.DataFrame(centers, columns=RATIO_COLUMNS))
    _, distances = pairwise_distances_argmin_min(train_scaled, center_scaled)
    coverage = pd.DataFrame(
        [
            {
                "K": background_size,
                "distance_mean": float(distances.mean()),
                "distance_median": float(np.median(distances)),
                "distance_p95": float(np.quantile(distances, 0.95)),
                "distance_max": float(distances.max()),
                "weight_sum": float(weights.sum()),
                "composition_sum_max_abs_error": float(composition_error.max()),
                "minimum_center_ratio": float(centers.min()),
            }
        ]
    )
    return center_frame, pd.DataFrame(moment_rows), coverage


def build_target_results(
    metals: pd.DataFrame,
    shap_values: np.ndarray,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    target_results = {}
    importance_columns = {}
    x_values = metals[RATIO_COLUMNS].reset_index(drop=True)
    for target in TARGETS:
        slug = target["slug"]
        output_index = int(target["index"])
        target_shap = shap_values[:, :, output_index]
        importance = importance_frame(target_shap, RATIO_COLUMNS)
        shap_frame = pd.DataFrame(target_shap, columns=RATIO_COLUMNS)
        average = pd.DataFrame(
            [importance["mean_abs_shap"].to_numpy()],
            columns=importance["feature"],
        )
        percent = pd.DataFrame(
            [importance["fraction"].to_numpy()],
            columns=importance["feature"],
        )
        target_results[slug] = {
            "x_val": x_values.copy(),
            "shap": shap_frame,
            "average": average,
            "percent": percent,
            "importance": importance,
        }
        importance_columns[target["physical_column"]] = pd.Series(
            np.abs(target_shap).mean(axis=0),
            index=RATIO_COLUMNS,
        )

    importance_matrix = pd.DataFrame(importance_columns).loc[
        RATIO_COLUMNS, DESCRIPTOR_PHYSICAL_COLUMNS
    ]
    direction_correlation_matrix, signed_effect_matrix = build_direction_matrices(
        x_values,
        target_results,
        importance_matrix,
    )
    return (
        target_results,
        importance_matrix,
        direction_correlation_matrix,
        signed_effect_matrix,
    )


def evaluate_background_stability(
    global_results: dict,
    descriptor_output_std: np.ndarray,
) -> tuple[int, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    importance_rows = []
    for background_size, result in global_results.items():
        importance_matrix = result["importance_matrix"]
        correlation_matrix = result["direction_correlation_matrix"]
        for output_index, descriptor in enumerate(DESCRIPTOR_PHYSICAL_COLUMNS):
            values = importance_matrix[descriptor]
            percentages = values / values.sum() * 100
            ranks = values.rank(method="min", ascending=False).astype(int)
            for metal in RATIO_COLUMNS:
                importance_rows.append(
                    {
                        "K": background_size,
                        "descriptor": descriptor,
                        "metal": metal,
                        "mean_abs_shap": float(values.loc[metal]),
                        "percent": float(percentages.loc[metal]),
                        "rank": int(ranks.loc[metal]),
                        "direction_correlation": float(correlation_matrix.loc[metal, descriptor]),
                        "signed_mean_abs_shap": float(
                            np.sign(correlation_matrix.loc[metal, descriptor])
                            * values.loc[metal]
                        ),
                        "expected_value": float(result["expected_value"][output_index]),
                    }
                )
    importance_long = pd.DataFrame(importance_rows)

    reference = global_results[GLOBAL_REFERENCE_BACKGROUND_SIZE]
    stability_rows = []
    candidate_pass = {}
    for background_size in GLOBAL_BACKGROUND_SIZES:
        descriptor_passes = []
        for output_index, descriptor in enumerate(DESCRIPTOR_PHYSICAL_COLUMNS):
            candidate_values = global_results[background_size]["importance_matrix"][descriptor]
            reference_values = reference["importance_matrix"][descriptor]
            candidate_percent = candidate_values / candidate_values.sum() * 100
            reference_percent = reference_values / reference_values.sum() * 100
            rank_correlation = float(
                candidate_percent.corr(reference_percent, method="spearman")
            )
            top_feature_same = bool(
                candidate_values.idxmax() == reference_values.idxmax()
            )
            max_percent_difference = float(
                np.max(np.abs(candidate_percent - reference_percent))
            )
            candidate_direction = global_results[background_size][
                "direction_correlation_matrix"
            ][descriptor]
            reference_direction = reference["direction_correlation_matrix"][descriptor]
            major_mask = reference_percent >= GLOBAL_MAJOR_IMPORTANCE_PERCENT
            direction_stable = bool(
                np.array_equal(
                    np.sign(candidate_direction.loc[major_mask]),
                    np.sign(reference_direction.loc[major_mask]),
                )
            )
            baseline_std_fraction = float(
                abs(
                    global_results[background_size]["expected_value"][output_index]
                    - reference["expected_value"][output_index]
                )
                / max(float(descriptor_output_std[output_index]), 1e-12)
            )
            passed = bool(
                rank_correlation >= GLOBAL_RANK_CORRELATION_THRESHOLD
                and top_feature_same
                and max_percent_difference <= GLOBAL_IMPORTANCE_PERCENT_TOLERANCE
                and direction_stable
                and baseline_std_fraction <= GLOBAL_BASELINE_STD_TOLERANCE
            )
            if background_size == GLOBAL_REFERENCE_BACKGROUND_SIZE:
                passed = True
            descriptor_passes.append(passed)
            stability_rows.append(
                {
                    "K": background_size,
                    "reference_K": GLOBAL_REFERENCE_BACKGROUND_SIZE,
                    "descriptor": descriptor,
                    "spearman_rank_correlation": rank_correlation,
                    "top_feature_same": top_feature_same,
                    "max_importance_difference_percentage_points": max_percent_difference,
                    "major_direction_stable": direction_stable,
                    "baseline_difference_fraction_of_output_std": baseline_std_fraction,
                    "passed": passed,
                }
            )
        candidate_pass[background_size] = bool(all(descriptor_passes))

    selected_size = GLOBAL_REFERENCE_BACKGROUND_SIZE
    for background_size in GLOBAL_BACKGROUND_SIZES:
        if candidate_pass[background_size]:
            selected_size = background_size
            break
    selection = pd.DataFrame(
        [
            {
                "selected_K": selected_size,
                "reference_K": GLOBAL_REFERENCE_BACKGROUND_SIZE,
                "tested_K": ", ".join(map(str, GLOBAL_BACKGROUND_SIZES)),
                "selection_rule": (
                    "Smallest K passing all descriptor-level rank, top-feature, "
                    "importance, direction, and baseline criteria"
                ),
                "selected_K_passed": candidate_pass[selected_size],
            }
        ]
    )
    return selected_size, importance_long, pd.DataFrame(stability_rows), selection


def save_background_stability_plot(importance_long: pd.DataFrame, path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    for output_index, (ax, descriptor) in enumerate(
        zip(axes.flat, DESCRIPTOR_PHYSICAL_COLUMNS)
    ):
        subset = importance_long.loc[importance_long["descriptor"] == descriptor]
        for metal in RATIO_COLUMNS:
            metal_values = subset.loc[subset["metal"] == metal].sort_values("K")
            ax.plot(
                metal_values["K"],
                metal_values["percent"],
                marker="o",
                linewidth=1.8,
                label=metal,
            )
        ax.set_title(DESCRIPTOR_PLOT_LABELS[output_index], fontsize=15)
        ax.set_ylabel("Mean |SHAP Value| (%)")
        ax.set_xticks(list(GLOBAL_BACKGROUND_SIZES))
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.8)
    axes[1, 0].set_xlabel("Background size K")
    axes[1, 1].set_xlabel("Background size K")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", frameon=False)
    fig.tight_layout(rect=[0, 0, 0.93, 1])
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def metadata_frame(metadata: dict) -> pd.DataFrame:
    return pd.DataFrame([metadata])


def run_descriptor_to_e10(
    train_metals: pd.DataFrame,
    test_metals: pd.DataFrame,
    train_descriptors: pd.DataFrame,
    test_descriptors: pd.DataFrame,
    assets: dict,
    metadata: dict,
):
    def model_fn(descriptor_values):
        return predict_e10_from_descriptors(descriptor_values, assets)

    background_values = train_descriptors[DESCRIPTOR_PHYSICAL_COLUMNS].to_numpy(dtype=float)
    x_values = test_descriptors[DESCRIPTOR_PHYSICAL_COLUMNS].to_numpy(dtype=float)
    shap_cube, expected_value = compute_kernel_shap(
        model_fn,
        background_values,
        x_values,
        n_outputs=1,
    )
    additivity = validate_additivity(
        model_fn,
        x_values,
        shap_cube,
        expected_value,
        ["E10"],
    )
    shap_values = shap_cube[:, :, 0]
    importance = importance_frame(shap_values, DESCRIPTOR_DISPLAY_COLUMNS)
    plot_importance = importance.copy()
    plot_importance["feature"] = plot_importance["feature"].map(DESCRIPTOR_DISPLAY_TO_PLOT)

    background_display = train_descriptors.copy()
    background_display.columns = DESCRIPTOR_DISPLAY_COLUMNS
    x_val = test_descriptors.copy()
    x_val.columns = DESCRIPTOR_DISPLAY_COLUMNS
    shap_frame = pd.DataFrame(shap_values, columns=DESCRIPTOR_DISPLAY_COLUMNS)
    average = pd.DataFrame(
        [importance["mean_abs_shap"].to_numpy()],
        columns=importance["feature"],
    )
    percent = pd.DataFrame(
        [importance["fraction"].to_numpy()],
        columns=importance["feature"],
    )
    predictions = pd.DataFrame(
        {
            "source_index": test_metals["source_index"].to_numpy(dtype=int),
            "predicted_E10": model_fn(x_values).reshape(-1),
            "reconstructed_E10": expected_value[0] + shap_values.sum(axis=1),
        }
    )
    stage_metadata = {
        **metadata,
        "analysis": "descriptor_to_E10",
        "background_source": "85 experimental training descriptor predictions",
        "explanation_source": "22 experimental test descriptor predictions",
        "background_count": int(len(background_values)),
        "explanation_count": int(len(x_values)),
        "kernel_nsamples": exact_kernel_nsamples(len(DESCRIPTOR_PHYSICAL_COLUMNS)),
        "target_unit": "physical E10 unit",
    }

    with pd.ExcelWriter(RESULT_DIR / "descriptor_to_E10_shap_data.xlsx") as writer:
        background_display.to_excel(writer, sheet_name="x_background", index=False)
        x_val.to_excel(writer, sheet_name="x_val", index=False)
        shap_frame.to_excel(writer, sheet_name="shap", index=False)
        average.to_excel(writer, sheet_name="average", index=False)
        percent.to_excel(writer, sheet_name="percent", index=False)
        importance.to_excel(writer, sheet_name="importance", index=False)
        train_descriptors.to_excel(
            writer,
            sheet_name="background_physical_raw",
            index=False,
        )
        test_descriptors.to_excel(
            writer,
            sheet_name="x_val_physical_raw",
            index=False,
        )
        train_metals.to_excel(writer, sheet_name="background_metals", index=False)
        test_metals.to_excel(writer, sheet_name="metals", index=False)
        predictions.to_excel(writer, sheet_name="predictions", index=False)
        additivity.to_excel(writer, sheet_name="additivity", index=False)
        metadata_frame(stage_metadata).to_excel(
            writer,
            sheet_name="metadata",
            index=False,
        )

    save_summary_plot(
        shap_values,
        x_values,
        DESCRIPTOR_PLOT_LABELS,
        RESULT_DIR / "descriptor_to_E10_shap_summary_plot.png",
    )
    save_combined_importance_plot(
        plot_importance,
        "Mean |SHAP Value|",
        RESULT_DIR / "descriptor_to_E10_shap_all.png",
    )
def build_direction_matrices(metals: pd.DataFrame, target_results: dict, importance_matrix: pd.DataFrame):
    corr_columns = {}
    signed_columns = {}
    for target in TARGETS:
        slug = target["slug"]
        descriptor = target["physical_column"]
        shap_values = target_results[slug]["shap"]
        corr = metals[RATIO_COLUMNS].corrwith(shap_values[RATIO_COLUMNS]).replace([np.inf, -np.inf], np.nan).fillna(0.0)
        signed_effect = np.sign(corr) * importance_matrix[descriptor]
        corr_columns[descriptor] = corr
        signed_columns[descriptor] = signed_effect
    direction_correlation_matrix = pd.DataFrame(corr_columns).loc[RATIO_COLUMNS, DESCRIPTOR_PHYSICAL_COLUMNS]
    signed_effect_matrix = pd.DataFrame(signed_columns).loc[RATIO_COLUMNS, DESCRIPTOR_PHYSICAL_COLUMNS]
    return direction_correlation_matrix, signed_effect_matrix


def save_metal_target_plots(
    target_results: dict,
    metals: pd.DataFrame,
    signed_effect_matrix: pd.DataFrame,
    prefix: str = "",
):
    x_values = metals[RATIO_COLUMNS].to_numpy(dtype=float)
    for target in TARGETS:
        slug = target["slug"]
        output_index = int(target["index"])
        target_label = DESCRIPTOR_PLOT_LABELS[output_index]
        save_summary_plot(
            target_results[slug]["shap"].to_numpy(dtype=float),
            x_values,
            RATIO_COLUMNS,
            RESULT_DIR / f"{prefix}metal_to_{slug}_shap_summary_plot.png",
        )
        save_combined_importance_plot(
            target_results[slug]["importance"],
            f"Mean |SHAP Value| for {target_label}",
            RESULT_DIR / f"{prefix}metal_to_{slug}_shap_signed_all.png",
            signed_effect=signed_effect_matrix[target["physical_column"]],
        )


def run_metal_to_descriptor_bridge(
    train_metals: pd.DataFrame,
    test_metals: pd.DataFrame,
    assets: dict,
    metadata: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    def model_fn(metal_values):
        return predict_descriptors_physical(metal_values, assets)

    background_values = train_metals[RATIO_COLUMNS].to_numpy(dtype=float)
    x_values = test_metals[RATIO_COLUMNS].to_numpy(dtype=float)
    train_descriptors = pd.DataFrame(
        model_fn(background_values),
        columns=DESCRIPTOR_PHYSICAL_COLUMNS,
    )
    test_descriptors = pd.DataFrame(
        model_fn(x_values),
        columns=DESCRIPTOR_PHYSICAL_COLUMNS,
    )
    shap_values, expected_value = compute_kernel_shap(
        model_fn,
        background_values,
        x_values,
        n_outputs=len(TARGETS),
    )
    additivity = validate_additivity(
        model_fn,
        x_values,
        shap_values,
        expected_value,
        DESCRIPTOR_PHYSICAL_COLUMNS,
    )
    (
        target_results,
        importance_matrix,
        direction_correlation_matrix,
        signed_effect_matrix,
    ) = build_target_results(test_metals, shap_values)

    save_metal_target_plots(
        target_results,
        test_metals,
        signed_effect_matrix,
    )
    save_descriptor_heatmaps(importance_matrix, signed_effect_matrix)

    expected_frame = pd.DataFrame(
        {
            "descriptor": DESCRIPTOR_PHYSICAL_COLUMNS,
            "expected_value": expected_value,
        }
    )
    stage_metadata = {
        **metadata,
        "analysis": "metal_to_descriptor_bridge",
        "background_source": "85 experimental training compositions",
        "explanation_source": "22 experimental test compositions",
        "background_count": int(len(background_values)),
        "explanation_count": int(len(x_values)),
        "kernel_nsamples": exact_kernel_nsamples(len(RATIO_COLUMNS)),
    }
    with pd.ExcelWriter(RESULT_DIR / "metal_to_descriptor_shap_data.xlsx") as writer:
        train_metals.to_excel(writer, sheet_name="background_metals", index=False)
        test_metals.to_excel(writer, sheet_name="metals", index=False)
        train_descriptors.to_excel(
            writer,
            sheet_name="background_descriptors",
            index=False,
        )
        test_descriptors.to_excel(
            writer,
            sheet_name="descriptor_predictions",
            index=False,
        )
        importance_matrix.to_excel(writer, sheet_name="importance_matrix")
        direction_correlation_matrix.to_excel(
            writer,
            sheet_name="direction_correlation_matrix",
        )
        signed_effect_matrix.to_excel(writer, sheet_name="signed_effect_matrix")
        expected_frame.to_excel(writer, sheet_name="expected_values", index=False)
        additivity.to_excel(writer, sheet_name="additivity", index=False)
        metadata_frame(stage_metadata).to_excel(
            writer,
            sheet_name="metadata",
            index=False,
        )
        for target in TARGETS:
            slug = target["slug"]
            result = target_results[slug]
            result["x_val"].to_excel(writer, sheet_name=f"x_val_{slug}", index=False)
            result["shap"].to_excel(writer, sheet_name=f"shap_{slug}", index=False)
            result["average"].to_excel(writer, sheet_name=f"average_{slug}", index=False)
            result["percent"].to_excel(writer, sheet_name=f"percent_{slug}", index=False)
            result["importance"].to_excel(
                writer,
                sheet_name=f"importance_{slug}",
                index=False,
            )
    return train_descriptors, test_descriptors
def run_global_metal_to_descriptor(
    theory_train: pd.DataFrame,
    theory_test: pd.DataFrame,
    assets: dict,
    metadata: dict,
) -> int:
    def model_fn(metal_values):
        return predict_descriptors_physical(metal_values, assets)

    train_values = theory_train[RATIO_COLUMNS]
    x_values = theory_test[RATIO_COLUMNS].to_numpy(dtype=float)
    train_predictions = model_fn(train_values.to_numpy(dtype=float))
    test_predictions = model_fn(x_values)
    descriptor_output_std = train_predictions.std(axis=0)

    global_results = {}
    center_frames = {}
    moment_frames = []
    coverage_frames = []
    for background_size in GLOBAL_BACKGROUND_SIZES:
        print(f"Running global metal-to-descriptor SHAP with K={background_size}...")
        background = shap.kmeans(
            train_values,
            background_size,
            round_values=False,
        )
        center_frame, moment_frame, coverage_frame = background_quality_frames(
            theory_train,
            background,
            assets["metal_scaler"],
            background_size,
        )
        shap_values, expected_value = compute_kernel_shap(
            model_fn,
            background,
            x_values,
            n_outputs=len(TARGETS),
        )
        additivity = validate_additivity(
            model_fn,
            x_values,
            shap_values,
            expected_value,
            DESCRIPTOR_PHYSICAL_COLUMNS,
        )
        additivity.insert(0, "K", background_size)
        (
            target_results,
            importance_matrix,
            direction_correlation_matrix,
            signed_effect_matrix,
        ) = build_target_results(theory_test, shap_values)
        global_results[background_size] = {
            "shap_values": shap_values,
            "expected_value": expected_value,
            "target_results": target_results,
            "importance_matrix": importance_matrix,
            "direction_correlation_matrix": direction_correlation_matrix,
            "signed_effect_matrix": signed_effect_matrix,
            "additivity": additivity,
        }
        center_frames[background_size] = center_frame
        moment_frames.append(moment_frame)
        coverage_frames.append(coverage_frame)

    selected_size, importance_long, stability, selection = (
        evaluate_background_stability(global_results, descriptor_output_std)
    )
    selected = global_results[selected_size]
    save_metal_target_plots(
        selected["target_results"],
        theory_test,
        selected["signed_effect_matrix"],
        prefix="global_",
    )
    save_descriptor_heatmaps(
        selected["importance_matrix"],
        selected["signed_effect_matrix"],
        prefix="global_",
    )
    save_background_stability_plot(
        importance_long,
        RESULT_DIR / "global_background_stability.png",
    )

    test_frame = theory_test.copy()
    for output_index, descriptor in enumerate(DESCRIPTOR_PHYSICAL_COLUMNS):
        test_frame[f"predicted_{descriptor}"] = test_predictions[:, output_index]
    all_additivity = pd.concat(
        [global_results[size]["additivity"] for size in GLOBAL_BACKGROUND_SIZES],
        ignore_index=True,
    )
    expected_rows = []
    for background_size in GLOBAL_BACKGROUND_SIZES:
        for output_index, descriptor in enumerate(DESCRIPTOR_PHYSICAL_COLUMNS):
            expected_rows.append(
                {
                    "K": background_size,
                    "descriptor": descriptor,
                    "expected_value": float(
                        global_results[background_size]["expected_value"][output_index]
                    ),
                    "training_prediction_std": float(descriptor_output_std[output_index]),
                }
            )
    global_metadata = {
        **metadata,
        "analysis": "global_metal_to_descriptor",
        "background_method": "weighted shap.kmeans with round_values=False",
        "background_sizes_tested": ", ".join(map(str, GLOBAL_BACKGROUND_SIZES)),
        "selected_background_size": int(selected_size),
        "explanation_count": int(len(theory_test)),
        "kernel_nsamples": exact_kernel_nsamples(len(RATIO_COLUMNS)),
        "shap_chunk_size": SHAP_CHUNK_SIZE,
    }

    with pd.ExcelWriter(
        RESULT_DIR / "metal_to_descriptor_global_shap_data.xlsx"
    ) as writer:
        test_frame.to_excel(writer, sheet_name="x_test", index=False)
        pd.concat(moment_frames, ignore_index=True).to_excel(
            writer,
            sheet_name="background_moments",
            index=False,
        )
        pd.concat(coverage_frames, ignore_index=True).to_excel(
            writer,
            sheet_name="background_coverage",
            index=False,
        )
        importance_long.to_excel(
            writer,
            sheet_name="importance_stability",
            index=False,
        )
        stability.to_excel(writer, sheet_name="selection_tests", index=False)
        selection.to_excel(writer, sheet_name="selection", index=False)
        pd.DataFrame(expected_rows).to_excel(
            writer,
            sheet_name="expected_values",
            index=False,
        )
        all_additivity.to_excel(writer, sheet_name="additivity", index=False)
        metadata_frame(global_metadata).to_excel(
            writer,
            sheet_name="metadata",
            index=False,
        )
        for background_size in GLOBAL_BACKGROUND_SIZES:
            center_frames[background_size].to_excel(
                writer,
                sheet_name=f"background_K{background_size}",
                index=False,
            )
            result = global_results[background_size]
            for target in TARGETS:
                slug = target["slug"]
                result["target_results"][slug]["shap"].to_excel(
                    writer,
                    sheet_name=f"shap_K{background_size}_{slug}",
                    index=False,
                )
    print(f"Selected global background size: K={selected_size}")
    return selected_size


def main():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    set_plot_style()
    assets = load_assets()
    theory_train, theory_test, theory_metadata = load_theory_split(assets)
    experiment_train, experiment_test, experiment_metadata = (
        load_experiment_split()
    )
    metadata = {
        **theory_metadata,
        **experiment_metadata,
        "weight_dir": repository_relative(WEIGHT_DIR),
        "theory_xlsx": repository_relative(THEORY_XLSX),
        "experiment_xlsx": repository_relative(EXPERIMENT_XLSX),
        "ratio_columns": ", ".join(RATIO_COLUMNS),
        "descriptor_physical_columns": ", ".join(DESCRIPTOR_PHYSICAL_COLUMNS),
        "original_combined_model_shap_rerun": False,
    }

    selected_size = run_global_metal_to_descriptor(
        theory_train,
        theory_test,
        assets,
        metadata,
    )
    train_descriptors, test_descriptors = run_metal_to_descriptor_bridge(
        experiment_train,
        experiment_test,
        assets,
        {
            **metadata,
            "selected_global_background_size": selected_size,
        },
    )
    run_descriptor_to_e10(
        experiment_train,
        experiment_test,
        train_descriptors,
        test_descriptors,
        assets,
        {
            **metadata,
            "selected_global_background_size": selected_size,
        },
    )
    print(f"Saved SHAP outputs to {RESULT_DIR}")
if __name__ == "__main__":
    main()
