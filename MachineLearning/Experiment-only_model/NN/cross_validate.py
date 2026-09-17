from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

from nn_common import (
    DirectExperimentModel,
    get_run_dirs_by_id,
    load_xy,
    make_loader,
    plot_cv_pred_vs_true,
    save_error_by_tier_outputs,
    save_json,
    train_model,
)


RUN_ID = "experiment_only_reference"

CONFIG = {
    "task": "direct_x_to_z_10fold_cv",
    "cv_seed": 106,
    "n_splits": 10,
    "batch_size": 8,
    "epochs": 300,
    "lr": 0.001,
    "optimizer": "Adamax",
    "scheduler_factor": 0.5,
    "scheduler_min_lr": 1e-6,
    "scheduler_patience": 10,
    "dropout": 0.4,
    "weight_decay": 1e-5,
    "input_dim": 5,
    "hidden_units": 512,
    "output_dim": 1,
}


def main():
    started_at = datetime.now()
    base_dir = Path(__file__).resolve().parent
    run_name, model_dir, result_dir = get_run_dirs_by_id(base_dir, RUN_ID)
    x, y = load_xy("experiment.xlsx", "metals", "potential")

    fold_rows = []
    pred_rows = []
    kfold = KFold(
        n_splits=CONFIG["n_splits"],
        shuffle=True,
        random_state=CONFIG["cv_seed"],
    )
    for fold, (train_idx, test_idx) in enumerate(kfold.split(x), start=1):
        x_train, x_test = x.iloc[train_idx], x.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        norm_x = StandardScaler().fit(x_train)
        norm_y = StandardScaler().fit(y_train)
        x_train_scaled = norm_x.transform(x_train)
        x_test_scaled = norm_x.transform(x_test)
        y_train_scaled = norm_y.transform(y_train)
        y_test_scaled = norm_y.transform(y_test)

        model = DirectExperimentModel(
            CONFIG["input_dim"],
            CONFIG["hidden_units"],
            CONFIG["output_dim"],
            CONFIG["dropout"],
        )
        train_model(
            model,
            make_loader(x_train_scaled, y_train_scaled, CONFIG["batch_size"]),
            make_loader(x_test_scaled, y_test_scaled, CONFIG["batch_size"]),
            CONFIG,
        )
        model.eval()
        with torch.no_grad():
            test_pred_scaled = model(
                torch.tensor(x_test_scaled).to(torch.float32)
            ).detach().numpy()
        test_pred = norm_y.inverse_transform(test_pred_scaled)
        y_true = y_test.values.reshape(-1)
        y_pred = test_pred.reshape(-1)
        corr = float(np.corrcoef(y_true, y_pred)[0, 1])
        rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
        fold_rows.append({"fold": fold, "corr": corr, "rmse": rmse})
        for row_id, true_value, pred_value in zip(test_idx, y_true, y_pred):
            pred_rows.append(
                {
                    "fold": fold,
                    "row_index": int(row_id),
                    "y_true": float(true_value),
                    "y_pred": float(pred_value),
                    "error": float(pred_value - true_value),
                    "abs_error": float(abs(pred_value - true_value)),
                }
            )

    fold_df = pd.DataFrame(fold_rows)
    pred_df = pd.DataFrame(pred_rows).sort_values("row_index")
    with pd.ExcelWriter(result_dir / "cv_fold_metrics.xlsx") as writer:
        fold_df.to_excel(writer, sheet_name="fold_metrics", index=False)
        pred_df.to_excel(writer, sheet_name="predictions", index=False)

    pred_metric = plot_cv_pred_vs_true(
        pred_df["y_true"],
        pred_df["y_pred"],
        result_dir / "cv_pred_vs_true.png",
        title="Direct x to z 10-Fold CV Prediction",
    )
    thresholds = save_error_by_tier_outputs(
        pred_df["y_true"], pred_df["y_pred"], result_dir, prefix="cv"
    )
    summary = {
        "run_name": run_name,
        "script": "cross_validate.py",
        "task": CONFIG["task"],
        "started_at": started_at.isoformat(timespec="seconds"),
        "ended_at": datetime.now().isoformat(timespec="seconds"),
        "config": CONFIG,
        "fold_corr_mean": float(fold_df["corr"].mean()),
        "fold_corr_std": float(fold_df["corr"].std(ddof=0)),
        "fold_rmse_mean": float(fold_df["rmse"].mean()),
        "fold_rmse_std": float(fold_df["rmse"].std(ddof=0)),
        "overall_metric": pred_metric,
        "error_tier_thresholds": thresholds,
    }
    save_json(summary, result_dir / "cv_summary.json")
    save_json(summary, model_dir / "cv_summary.json")
    print(f"Saved CV outputs to {result_dir}")


if __name__ == "__main__":
    main()
