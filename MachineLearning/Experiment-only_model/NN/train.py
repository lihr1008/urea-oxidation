from datetime import datetime
from pathlib import Path

import torch
from nn_common import (
    DirectExperimentModel,
    append_run_metadata,
    get_run_dirs_by_id,
    load_xy,
    make_loader,
    metrics,
    plot_loss_curve,
    plot_single,
    predict_inverse,
    save_error_by_tier_outputs,
    save_pickle,
    save_prediction_excel,
    save_run_metadata,
    split_scale,
    train_model,
)

RUN_ID = "experiment_only_reference"

CONFIG = {
    "task": "direct_x_to_z",
    "seed": 106,  # random seed
    "test_size": 0.2,  # test set proportion
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
    x, y = load_xy("experiment.xlsx", "metals", "potential")
    data = split_scale(x, y, CONFIG["seed"], CONFIG["test_size"])

    train_loader = make_loader(
        data["x_train_scaled"], data["y_train_scaled"], CONFIG["batch_size"]
    )
    test_loader = make_loader(
        data["x_test_scaled"], data["y_test_scaled"], CONFIG["batch_size"]
    )
    model = DirectExperimentModel(
        CONFIG["input_dim"],
        CONFIG["hidden_units"],
        CONFIG["output_dim"],
        CONFIG["dropout"],
    )

    loss_history = train_model(model, train_loader, test_loader, CONFIG)
    train_pred = predict_inverse(model, data["x_train_scaled"], data["norm_y"])
    test_pred = predict_inverse(model, data["x_test_scaled"], data["norm_y"])
    metric = metrics(data["y_train"], data["y_test"], train_pred, test_pred)
    print(metric)

    ended_at = datetime.now()
    run_name, model_dir, result_dir = get_run_dirs_by_id(base_dir, RUN_ID)
    torch.save(model.state_dict(), model_dir / "exp_model_weights.pth")
    torch.save(model, model_dir / "exp_model.pth")
    save_pickle(data["norm_x"], model_dir / "norm_x.pckl")
    save_pickle(data["norm_y"], model_dir / "norm_y.pckl")
    save_run_metadata(
        model_dir,
        CONFIG,
        {
            "run_name": run_name,
            "script": "train.py",
            "task": CONFIG["task"],
            "started_at": started_at.isoformat(timespec="seconds"),
            "ended_at": ended_at.isoformat(timespec="seconds"),
            "data_file": "MachineLearning/data/experiment.xlsx",
            "x_shape": list(x.shape),
            "y_shape": list(y.shape),
            "input_dim": CONFIG["input_dim"],
            "output_dim": CONFIG["output_dim"],
            "metrics": metric,
            "loss_history": loss_history,
        },
    )

    plot_single(
        y,
        data["y_train"],
        data["y_test"],
        train_pred,
        test_pred,
        metric,
        result_dir / "train.png",
    )
    plot_loss_curve(loss_history, result_dir / "loss_curve.png")
    error_thresholds = save_error_by_tier_outputs(
        data["y_test"], test_pred, result_dir
    )
    append_run_metadata(
        model_dir,
        {
            "error_tier_thresholds": error_thresholds,
            "loss_curve": "loss_curve.png",
            "error_distribution_violin_box": "error_distribution_violin_box.png",
            "error_distribution_butterfly": "error_distribution_butterfly.png",
        },
    )
    save_prediction_excel(
        result_dir / "result.xlsx",
        data["y_train"],
        train_pred,
        data["y_test"],
        test_pred,
        ["potential"],
    )
    print(f"Saved model information to {model_dir}")
    print(f"Saved run results to {result_dir}")


if __name__ == "__main__":
    main()
