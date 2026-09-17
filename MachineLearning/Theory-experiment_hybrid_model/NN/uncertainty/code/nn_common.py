import json
import pickle
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import KFold, train_test_split
from sklearn.preprocessing import StandardScaler
from torch.nn import functional as F
from torch.utils.data import DataLoader, Dataset


CODE_DIR = Path(__file__).resolve().parent
ML_ROOT = Path(__file__).resolve().parents[4]
REPO_ROOT = ML_ROOT.parent
DATA_DIR = ML_ROOT / "data"
RATIO_COLUMNS = ["Mn", "Fe", "Co", "Ni", "Zn"]
DESCRIPTOR_COLUMNS = ["ade_urea", "ade_CO", "e_urea", "e_CO"]
if str(ML_ROOT) not in sys.path:
    sys.path.insert(0, str(ML_ROOT))

from composition_ilr_pca import (
    fit_reference_ilr_pca,
    save_ilr_pca_outputs,
    transform_ilr_pca,
)


class ArrayDataset(Dataset):
    def __init__(self, data, label):
        self.data = data
        self.label = label

    def __getitem__(self, index):
        return self.data[index], self.label[index]

    def __len__(self):
        return len(self.data)


class DirectExperimentModel(nn.Module):
    def __init__(self, input_dim=5, hidden_units=64, output_dim=1):
        super().__init__()
        self.layer1 = nn.Linear(input_dim, hidden_units)
        self.layer2 = nn.Linear(hidden_units, hidden_units)
        self.layer3 = nn.Linear(hidden_units, hidden_units)
        self.layer4 = nn.Linear(hidden_units, output_dim)

    def forward(self, x):
        x = F.relu(self.layer1(x))
        x = F.relu(self.layer2(x))
        x = F.relu(self.layer3(x))
        return self.layer4(x)


class DescriptorPreModel(nn.Module):
    def __init__(self, input_dim=5, hidden_units=512, output_dim=4):
        super().__init__()
        self.layer1 = nn.Linear(input_dim, hidden_units)
        self.layer2 = nn.Linear(hidden_units, hidden_units)
        self.layer3 = nn.Linear(hidden_units, output_dim)

    def forward(self, x):
        x = F.relu(self.layer1(x))
        x = F.relu(self.layer2(x))
        return self.layer3(x)


class DescriptorToExperimentModel(nn.Module):
    def __init__(self, input_dim=4, hidden_units=512, output_dim=1, dropout=0.4):
        super().__init__()
        self.layer1 = nn.Linear(input_dim, hidden_units)
        self.layer2 = nn.Linear(hidden_units, hidden_units)
        self.layer3 = nn.Linear(hidden_units, hidden_units)
        self.layer4 = nn.Linear(hidden_units, output_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        x = self.dropout(F.relu(self.layer1(x)))
        x = self.dropout(F.relu(self.layer2(x)))
        x = self.dropout(F.relu(self.layer3(x)))
        return self.layer4(x)


def load_xy(workbook, x_sheet, y_sheet):
    excel = DATA_DIR / workbook
    available_sheets = set(pd.ExcelFile(excel).sheet_names)
    if y_sheet not in available_sheets:
        target_alias = {"potential": "overpotential", "overpotential": "potential"}.get(y_sheet)
        if target_alias in available_sheets:
            y_sheet = target_alias
    return (
        pd.read_excel(excel, sheet_name=x_sheet),
        pd.read_excel(excel, sheet_name=y_sheet),
    )


def split_scale(x, y, seed, test_size=0.2):
    x_train, x_test, y_train, y_test = train_test_split(
        x, y, test_size=test_size, random_state=int(seed)
    )
    norm_x = StandardScaler().fit(x_train)
    norm_y = StandardScaler().fit(y_train)
    return {
        "x_train": x_train,
        "x_test": x_test,
        "y_train": y_train,
        "y_test": y_test,
        "norm_x": norm_x,
        "norm_y": norm_y,
        "x_train_scaled": norm_x.transform(x_train),
        "x_test_scaled": norm_x.transform(x_test),
        "y_train_scaled": norm_y.transform(y_train),
        "y_test_scaled": norm_y.transform(y_test),
    }


def make_loader(x_data, y_data, batch_size):
    dataset = ArrayDataset(x_data, y_data)
    return DataLoader(dataset, batch_size=int(batch_size), shuffle=False)


def make_optimizer(model, config):
    name = config.get("optimizer", "Adam")
    lr = config.get("lr", 0.001)
    weight_decay = config.get("weight_decay", 0.0)
    if name == "Adamax":
        return optim.Adamax(model.parameters(), lr=lr, weight_decay=weight_decay)
    if name == "Adam":
        return optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    raise ValueError(f"Unsupported optimizer: {name}")


def make_scheduler(optimizer, config):
    return torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        "min",
        factor=config.get("scheduler_factor", 0.5),
        min_lr=config.get("scheduler_min_lr", 1e-6),
        patience=config.get("scheduler_patience", 10),
    )


def train_model(model, train_loader, test_loader, config):
    criterion = nn.MSELoss()
    optimizer = make_optimizer(model, config)
    scheduler = make_scheduler(optimizer, config)
    losses = []

    for epoch in range(int(config["epochs"])):
        model.train()
        train_losses = []
        for data, label in train_loader:
            data = data.to(torch.float32)
            label = label.to(torch.float32)
            logits = model(data)
            loss = criterion(logits, label)
            train_losses.append(loss.detach().item())
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        model.eval()
        test_losses = []
        with torch.no_grad():
            for data, label in test_loader:
                data = data.to(torch.float32)
                label = label.to(torch.float32)
                logits = model(data)
                loss = criterion(logits, label)
                test_losses.append(loss.detach().item())

        train_loss = float(np.average(train_losses))
        test_loss = float(np.average(test_losses))
        scheduler.step(train_loss)
        losses.append({"epoch": epoch, "train_loss": train_loss, "test_loss": test_loss})
        print(f"Train Epoch : {epoch}\tLoss:{train_loss}")
        print(f"\nTest epoch:{epoch}, Average loss: {test_loss}\n")

    return losses


def predict_inverse(model, x_scaled, norm_y):
    model.eval()
    with torch.no_grad():
        x_tensor = torch.tensor(x_scaled).to(torch.float32)
        pred = model(x_tensor).detach().numpy()
    return norm_y.inverse_transform(pred)


def metrics(y_train, y_test, train_pred, test_pred):
    train_pred = np.asarray(train_pred)
    test_pred = np.asarray(test_pred)
    y_train = np.asarray(y_train)
    y_test = np.asarray(y_test)
    train_corr = [
        float(np.corrcoef(y_train[:, i], train_pred[:, i])[0, 1])
        for i in range(train_pred.shape[1])
    ]
    test_corr = [
        float(np.corrcoef(y_test[:, i], test_pred[:, i])[0, 1])
        for i in range(test_pred.shape[1])
    ]
    train_rmse = [
        float(np.sqrt(np.mean((y_train[:, i] - train_pred[:, i]) ** 2)))
        for i in range(train_pred.shape[1])
    ]
    test_rmse = [
        float(np.sqrt(np.mean((y_test[:, i] - test_pred[:, i]) ** 2)))
        for i in range(test_pred.shape[1])
    ]
    return {
        "train_corr": train_corr,
        "test_corr": test_corr,
        "train_rmse": train_rmse,
        "test_rmse": test_rmse,
    }


def plot_single(y, y_train, y_test, train_pred, test_pred, metric, path, figsize=(9, 4)):
    fig, axes = plt.subplots(1, 2, figsize=figsize, tight_layout=True)
    y_min = min(y.values.min(), train_pred.min(), test_pred.min())
    y_max = max(y.values.max(), train_pred.max(), test_pred.max())
    lim = [y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05]

    for axis in axes:
        axis.set_aspect("equal")
        axis.plot(lim, lim, lw=1)
        axis.set_xlim(lim)
        axis.set_ylim(lim)

    axes[0].scatter(y_train.iloc[:, 0], train_pred[:, 0], s=5)
    axes[1].scatter(y_test.iloc[:, 0], test_pred[:, 0], s=5)
    axes[0].text(
        y_max,
        y_min,
        f"r={metric['train_corr'][0]:.4f}\nRMSE={metric['train_rmse'][0]:.4f}",
        horizontalalignment="right",
        fontsize=12,
    )
    axes[1].text(
        y_max,
        y_min,
        f"r={metric['test_corr'][0]:.4f}\nRMSE={metric['test_rmse'][0]:.4f}",
        horizontalalignment="right",
        fontsize=12,
    )
    axes[0].set_ylabel(y.columns[0], fontsize=16)
    fig.savefig(path, dpi=300, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)


def plot_multi(y, y_train, y_test, train_pred, test_pred, metric, path):
    fig, axes = plt.subplots(y.shape[1], 2, figsize=(9, 12), tight_layout=True)
    axes[0, 0].set_title("Train", fontsize=16)
    axes[0, 1].set_title("Test", fontsize=16)
    for i in range(y.shape[1]):
        y_min = min(y.iloc[:, i].min(), train_pred[:, i].min(), test_pred[:, i].min())
        y_max = max(y.iloc[:, i].max(), train_pred[:, i].max(), test_pred[:, i].max())
        lim = [y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05]
        for j in [0, 1]:
            axes[i, j].set_aspect("equal")
            axes[i, j].plot(lim, lim, lw=1)
            axes[i, j].set_xlim(lim)
            axes[i, j].set_ylim(lim)
        axes[i, 0].scatter(y_train.iloc[:, i], train_pred[:, i], s=5)
        axes[i, 1].scatter(y_test.iloc[:, i], test_pred[:, i], s=5)
        axes[i, 0].text(
            y_max,
            y_min,
            f"r={metric['train_corr'][i]:.4f}\nRMSE={metric['train_rmse'][i]:.4f}",
            horizontalalignment="right",
            fontsize=12,
        )
        axes[i, 1].text(
            y_max,
            y_min,
            f"r={metric['test_corr'][i]:.4f}\nRMSE={metric['test_rmse'][i]:.4f}",
            horizontalalignment="right",
            fontsize=12,
        )
        axes[i, 0].set_ylabel(y.columns[i], fontsize=16)
    fig.savefig(path, dpi=300, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)


def save_prediction_excel(path, y_train, train_pred, y_test, test_pred, columns):
    with pd.ExcelWriter(path) as writer:
        y_train.to_excel(writer, sheet_name="y_train", index=False)
        pd.DataFrame(train_pred, columns=columns).to_excel(
            writer, sheet_name="train_pred", index=False
        )
        y_test.to_excel(writer, sheet_name="y_test", index=False)
        pd.DataFrame(test_pred, columns=columns).to_excel(
            writer, sheet_name="test_pred", index=False
        )


def save_pickle(obj, path):
    with open(path, "wb") as f:
        pickle.dump(obj, f)


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def save_json(data, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def save_run_metadata(model_dir, config, run_info):
    save_json(config, model_dir / "config.json")
    save_json(run_info, model_dir / "run_info.json")


def append_run_metadata(model_dir, updates):
    path = Path(model_dir) / "run_info.json"
    if path.exists():
        with open(path, encoding="utf-8") as f:
            run_info = json.load(f)
    else:
        run_info = {}
    run_info.update(updates)
    save_json(run_info, path)


def read_choose_proportion():
    rows = []
    with open(DATA_DIR / "choose_proportion.txt") as f:
        next(f)
        for line in f:
            if line.strip():
                rows.append([float(value) / 100 for value in line.split()])
    return np.asarray(rows)


def save_space_prediction(path, x_val, values):
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame(x_val, columns=RATIO_COLUMNS).to_excel(
            writer, sheet_name="ratio", index=False
        )
        pd.DataFrame(values, columns=["potential"]).to_excel(
            writer, sheet_name="value", index=False
        )


def create_analysis_from_space_prediction(space_path, analysis_path):
    space_path = Path(space_path)
    analysis_path = Path(analysis_path)
    if not space_path.exists():
        raise FileNotFoundError(
            f"{space_path} does not exist. Run grid_search.py before dimension_reduction.py."
        )

    ratio = pd.read_excel(space_path, sheet_name="ratio")
    value = pd.read_excel(space_path, sheet_name="value")
    if ratio.shape[1] < 5:
        raise ValueError(f"{space_path} ratio sheet must contain Mn/Fe/Co/Ni/Zn columns.")
    if value.shape[1] < 1:
        raise ValueError(f"{space_path} value sheet must contain a predicted potential column.")

    ratio = ratio.iloc[:, :5].copy()
    ratio.columns = RATIO_COLUMNS
    potential = value.iloc[:, 0].copy()
    analysis = ratio.copy()
    analysis["potential"] = potential
    analysis = analysis.sort_values("potential", ascending=True).reset_index(drop=True)
    analysis.to_excel(analysis_path, index=False)
    return analysis


def run_ilr_pca(input_path, output_dir):
    data = pd.read_excel(input_path)
    validate_analysis_frame(data, input_path)
    return save_ilr_pca_outputs(
        data,
        output_dir,
        DATA_DIR / "experiment.xlsx",
        target_column=data.columns[-1],
        target_label="potential",
    )


def plot_loss_curve(loss_history, path, title="Training and Test Loss"):
    if not loss_history:
        return
    history = pd.DataFrame(loss_history)
    fig, ax = plt.subplots(figsize=(8, 5), tight_layout=True)
    ax.plot(history["epoch"], history["train_loss"], label="Train loss", lw=2)
    ax.plot(history["epoch"], history["test_loss"], label="Test loss", lw=2)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("MSE loss")
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=300, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)


def error_by_tier(y_true, y_pred):
    y_true = np.asarray(y_true).reshape(-1)
    y_pred = np.asarray(y_pred).reshape(-1)
    q1, q2 = np.quantile(y_true, [1 / 3, 2 / 3])
    tiers = np.where(y_true <= q1, "low", np.where(y_true <= q2, "medium", "high"))
    return pd.DataFrame(
        {
            "y_true": y_true,
            "y_pred": y_pred,
            "error": y_pred - y_true,
            "abs_error": np.abs(y_pred - y_true),
            "tier": tiers,
        }
    ), {"low_medium": float(q1), "medium_high": float(q2)}


def save_error_by_tier_outputs(y_true, y_pred, output_dir, prefix=""):
    output_dir = Path(output_dir)
    data, thresholds = error_by_tier(y_true, y_pred)
    stem = f"{prefix}_" if prefix else ""
    data.to_excel(output_dir / f"{stem}error_by_tier.xlsx", index=False)
    plot_error_distribution_violin_box(
        data, output_dir / f"{stem}error_distribution_violin_box.png"
    )
    plot_error_distribution_butterfly(
        data, output_dir / f"{stem}error_distribution_butterfly.png"
    )
    return thresholds


def plot_error_distribution_violin_box(data, path):
    tiers = ["low", "medium", "high"]
    values = [data.loc[data["tier"] == tier, "error"].values for tier in tiers]
    fig, ax = plt.subplots(figsize=(8, 5), tight_layout=True)
    ax.violinplot(values, showmeans=False, showmedians=False, showextrema=False)
    ax.boxplot(values, widths=0.18, patch_artist=True, showfliers=True)
    ax.axhline(0, color="black", lw=1, alpha=0.6)
    ax.set_xticks(range(1, len(tiers) + 1))
    ax.set_xticklabels(tiers)
    ax.set_xlabel("Experimental performance tier")
    ax.set_ylabel("Prediction error (predicted - true)")
    ax.set_title("Error Distribution by Tier (Violin + Box Plot)")
    fig.savefig(path, dpi=300, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)


def plot_error_distribution_butterfly(data, path):
    tiers = ["low", "medium", "high"]
    fig, ax = plt.subplots(figsize=(8, 5), tight_layout=True)
    for idx, tier in enumerate(tiers):
        values = data.loc[data["tier"] == tier, "error"].values
        if len(values) == 0:
            continue
        counts, bins = np.histogram(values, bins=12)
        centers = (bins[:-1] + bins[1:]) / 2
        widths = counts / counts.max() * 0.35 if counts.max() else counts
        ax.barh(centers, widths, height=np.diff(bins), left=idx, alpha=0.55)
        ax.barh(centers, -widths, height=np.diff(bins), left=idx, alpha=0.55)
    ax.axhline(0, color="black", lw=1, alpha=0.6)
    ax.set_xticks(range(len(tiers)))
    ax.set_xticklabels(tiers)
    ax.set_xlabel("Experimental performance tier")
    ax.set_ylabel("Prediction error (predicted - true)")
    ax.set_title("Error Distribution by Tier (Butterfly Plot)")
    fig.savefig(path, dpi=300, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)


def plot_cv_pred_vs_true(y_true, y_pred, path, title="10-Fold CV Prediction"):
    y_true = np.asarray(y_true).reshape(-1)
    y_pred = np.asarray(y_pred).reshape(-1)
    y_min = min(y_true.min(), y_pred.min())
    y_max = max(y_true.max(), y_pred.max())
    lim = [y_min - (y_max - y_min) * 0.05, y_max + (y_max - y_min) * 0.05]
    corr = float(np.corrcoef(y_true, y_pred)[0, 1])
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    fig, ax = plt.subplots(figsize=(6, 6), tight_layout=True)
    ax.scatter(y_true, y_pred, s=18, alpha=0.75)
    ax.plot(lim, lim, color="red", lw=1)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Actual")
    ax.set_ylabel("Predicted")
    ax.set_title(title)
    ax.text(lim[1], lim[0], f"r={corr:.4f}\nRMSE={rmse:.4f}", ha="right", va="bottom")
    fig.savefig(path, dpi=300, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)
    return {"corr": corr, "rmse": rmse}


def composition_space_ilr_pca(output_dir, analysis_path=None):
    output_dir = Path(output_dir)
    grid = pd.DataFrame(read_choose_proportion(), columns=RATIO_COLUMNS)
    experiment_x, _ = load_xy("experiment.xlsx", "metals", "potential")
    nonzero_counts = (experiment_x[RATIO_COLUMNS] > 0).sum(axis=1)
    single_element = experiment_x.loc[nonzero_counts == 1, RATIO_COLUMNS].copy()
    initial = experiment_x.loc[nonzero_counts != 1, RATIO_COLUMNS].copy()

    top10 = pd.DataFrame(columns=RATIO_COLUMNS)
    top10_values = pd.Series(dtype=float)
    if analysis_path is not None and Path(analysis_path).exists():
        analysis = pd.read_excel(analysis_path)
        validate_analysis_frame(analysis, analysis_path)
        ratio_cols = list(analysis.columns[:5])
        value_col = analysis.columns[-1]
        top10_raw = analysis.sort_values(value_col, ascending=True).head(10)
        top10 = top10_raw.loc[:, ratio_cols].copy()
        top10.columns = RATIO_COLUMNS
        top10_values = top10_raw[value_col].reset_index(drop=True)

    pca = fit_reference_ilr_pca(DATA_DIR / "experiment.xlsx")

    def transform_frame(frame, label, values=None):
        if frame.empty:
            return pd.DataFrame(
                columns=RATIO_COLUMNS
                + ["ILR_PC1", "ILR_PC2", "source", "predicted_z"]
            )
        coords = transform_ilr_pca(pca, frame[RATIO_COLUMNS])
        out = frame[RATIO_COLUMNS].copy().reset_index(drop=True)
        out["ILR_PC1"] = coords[:, 0]
        out["ILR_PC2"] = coords[:, 1]
        out["source"] = label
        if values is not None and len(values) == len(out):
            out["predicted_z"] = list(values)
        else:
            out["predicted_z"] = np.nan
        return out

    grid_out = transform_frame(grid, "1% composition grid")
    initial_out = transform_frame(initial, "initial experimental points")
    single_out = transform_frame(single_element, "single-element points")
    top10_out = transform_frame(top10, "top-10 predicted candidates", top10_values)
    frames = [grid_out, initial_out, single_out]
    if not top10_out.empty:
        frames.append(top10_out)
    all_out = pd.concat(frames, ignore_index=True)
    all_out.to_excel(output_dir / "composition_space_ilr_pca.xlsx", index=False)

    fig, ax = plt.subplots(figsize=(9, 7), tight_layout=True)
    ax.scatter(
        grid_out["ILR_PC1"],
        grid_out["ILR_PC2"],
        s=1,
        color="lightgray",
        alpha=0.35,
        label="1% composition grid",
    )
    ax.scatter(
        initial_out["ILR_PC1"],
        initial_out["ILR_PC2"],
        s=32,
        color="#1f77b4",
        edgecolor="white",
        linewidth=0.4,
        label="initial experimental points",
    )
    ax.scatter(
        single_out["ILR_PC1"],
        single_out["ILR_PC2"],
        s=56,
        marker="s",
        color="#ff7f0e",
        edgecolor="black",
        linewidth=0.5,
        label="single-element points",
    )
    if not top10_out.empty:
        ax.scatter(
            top10_out["ILR_PC1"],
            top10_out["ILR_PC2"],
            s=72,
            marker="*",
            color="#d62728",
            edgecolor="black",
            linewidth=0.5,
            label="top-10 predicted candidates",
        )
    ax.set_xlabel("ILR-PC1")
    ax.set_ylabel("ILR-PC2")
    ax.set_title("Composition space ILR-PCA")
    ax.legend(markerscale=2)
    fig.savefig(
        output_dir / "composition_space_ilr_pca.png",
        dpi=300,
        bbox_inches="tight",
        metadata={"Software": None},
    )
    plt.close(fig)
    return {
        "grid_points": int(len(grid_out)),
        "initial_experimental_points": int(len(initial_out)),
        "single_element_points": int(len(single_out)),
        "top10_points": int(len(top10_out)),
    }


def validate_analysis_frame(data, source):
    if data.shape[1] < 6:
        raise ValueError(
            f"{source} must contain at least 6 columns: Mn, Fe, Co, Ni, Zn, potential."
        )
    first_five = list(data.columns[:5])
    if first_five != RATIO_COLUMNS:
        raise ValueError(
            f"{source} first five columns must be {RATIO_COLUMNS}, got {first_five}."
        )
