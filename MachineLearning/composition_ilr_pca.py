from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


RATIO_COLUMNS = ["Mn", "Fe", "Co", "Ni", "Zn"]
DEFAULT_ZERO_REPLACEMENT_DELTA = 0.00025
REFERENCE_COMPOSITION_COUNT = 102


def compositions_to_unit(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(RATIO_COLUMNS):
        raise ValueError(
            f"Compositions must be a two-dimensional array with {len(RATIO_COLUMNS)} columns."
        )
    if not np.isfinite(values).all():
        raise ValueError("Compositions must contain only finite values.")
    if np.any(values < 0):
        raise ValueError("Compositions cannot contain negative values.")
    row_sum = values.sum(axis=1, keepdims=True)
    if np.any(row_sum <= 0):
        raise ValueError("Composition rows must have positive sums before ILR transform.")
    return values / row_sum


def multiplicative_zero_replacement(values, delta=DEFAULT_ZERO_REPLACEMENT_DELTA):
    comp = compositions_to_unit(values)
    out = comp.copy()
    delta = float(delta)
    if not 0 < delta < 1:
        raise ValueError("The ILR zero-replacement delta must lie between 0 and 1.")

    for row_index in range(out.shape[0]):
        row = out[row_index].copy()
        zero_mask = row <= 0
        zero_count = int(zero_mask.sum())
        if zero_count == 0:
            continue
        replacement_total = zero_count * delta
        if replacement_total >= 1.0:
            raise ValueError("The ILR zero-replacement delta is too large.")
        nonzero_sum = float(row[~zero_mask].sum())
        row[zero_mask] = delta
        row[~zero_mask] *= (1.0 - replacement_total) / nonzero_sum
        out[row_index] = row
    return out


def ilr_basis(dim):
    basis = np.zeros((dim - 1, dim), dtype=float)
    for index in range(1, dim):
        basis[index - 1, :index] = 1.0 / np.sqrt(index * (index + 1))
        basis[index - 1, index] = -index / np.sqrt(index * (index + 1))
    return basis


def ilr_transform(values, delta=DEFAULT_ZERO_REPLACEMENT_DELTA):
    comp = multiplicative_zero_replacement(values, delta)
    return np.log(comp) @ ilr_basis(comp.shape[1]).T


def fit_reference_ilr_pca(
    experiment_path,
    delta=DEFAULT_ZERO_REPLACEMENT_DELTA,
):
    experiment_path = Path(experiment_path)
    experiment = pd.read_excel(experiment_path, sheet_name="metals")
    missing = [column for column in RATIO_COLUMNS if column not in experiment.columns]
    if missing:
        raise ValueError(
            f"{experiment_path} is missing composition columns: {', '.join(missing)}"
        )

    experiment_ratios = experiment[RATIO_COLUMNS].to_numpy(dtype=float)
    multicomponent_mask = (experiment_ratios > 0).sum(axis=1) > 1
    reference_count = int(multicomponent_mask.sum())
    if reference_count != REFERENCE_COMPOSITION_COUNT:
        raise ValueError(
            "The ILR-PCA reference must contain exactly "
            f"{REFERENCE_COMPOSITION_COUNT} multicomponent experimental compositions; "
            f"found {reference_count}."
        )

    experiment_ilr = ilr_transform(experiment_ratios, delta)
    pca = PCA(n_components=2, svd_solver="full")
    pca.fit(experiment_ilr[multicomponent_mask])

    # PCA signs are otherwise arbitrary. Match the convention used by the
    # reliability-map workflow by making the largest element loading positive.
    basis = ilr_basis(len(RATIO_COLUMNS))
    element_loadings = pca.components_ @ basis
    for component_index in range(pca.n_components_):
        pivot = int(np.argmax(np.abs(element_loadings[component_index])))
        if element_loadings[component_index, pivot] < 0:
            pca.components_[component_index] *= -1.0
    return pca


def transform_ilr_pca(pca, compositions, delta=DEFAULT_ZERO_REPLACEMENT_DELTA):
    if isinstance(compositions, pd.DataFrame):
        missing = [column for column in RATIO_COLUMNS if column not in compositions.columns]
        if missing:
            raise ValueError(
                f"Composition table is missing columns: {', '.join(missing)}"
            )
        values = compositions[RATIO_COLUMNS].to_numpy(dtype=float)
    else:
        values = np.asarray(compositions, dtype=float)
    return pca.transform(ilr_transform(values, delta))


def save_ilr_pca_outputs(
    data,
    output_dir,
    experiment_path,
    target_column=None,
    target_label="potential",
    delta=DEFAULT_ZERO_REPLACEMENT_DELTA,
):
    missing = [column for column in RATIO_COLUMNS if column not in data.columns]
    if missing:
        raise ValueError(f"Analysis table is missing columns: {', '.join(missing)}")
    if target_column is None:
        target_column = data.columns[-1]
    if target_column not in data.columns:
        raise ValueError(f"Analysis table is missing target column: {target_column}")

    pca = fit_reference_ilr_pca(experiment_path, delta)
    coordinates = transform_ilr_pca(pca, data[RATIO_COLUMNS], delta)
    target = data[target_column].to_numpy(dtype=float)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    result = pd.DataFrame(
        {
            "ILR_PC1": coordinates[:, 0],
            "ILR_PC2": coordinates[:, 1],
            target_label: target,
        }
    )
    result.to_csv(output_dir / "ilr_pca_results.csv", index=False)

    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["axes.unicode_minus"] = False
    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(
        coordinates[:, 0], coordinates[:, 1], c=target, cmap="jet"
    )
    colorbar = fig.colorbar(scatter, ax=ax)
    colorbar.set_label(target_label, fontsize=25)
    colorbar.ax.tick_params(labelsize=20)
    ax.set_title(
        "ILR-PCA fitted on 102 multicomponent experimental compositions",
        fontsize=19,
    )
    ax.set_xlabel("ILR-PC1", fontsize=25)
    ax.set_ylabel("ILR-PC2", fontsize=25)
    ax.tick_params(labelsize=20)
    fig.savefig(output_dir / "ilr_pca.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return result
