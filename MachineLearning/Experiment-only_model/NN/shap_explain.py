from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import shap
import torch
import matplotlib.pyplot as plt

from nn_common import (
    DirectExperimentModel,
    get_run_dirs_by_id,
    load_pickle,
    load_xy,
    split_scale,
)

RUN_ID = "experiment_only_reference"

def main(run_id=RUN_ID):
    base_dir = Path(__file__).resolve().parent
    run_name, model_dir, result_dir = get_run_dirs_by_id(base_dir, run_id)
    with open(model_dir / "config.json", encoding="utf-8") as f:
        config = json.load(f)

    model = DirectExperimentModel(
        int(config["input_dim"]),
        int(config["hidden_units"]),
        int(config["output_dim"]),
        float(config.get("dropout", 0.4)),
    )
    model.load_state_dict(torch.load(model_dir / "exp_model_weights.pth", map_location="cpu"))
    model.eval()
    norm_x = load_pickle(model_dir / "norm_x.pckl")

    x, y = load_xy("experiment.xlsx", "metals", "potential")
    data = split_scale(x, y, int(config["seed"]), config.get("test_size", 0.2))

    def combined_model(x_np):
        x_norm = norm_x.transform(x_np)
        x_tensor = torch.tensor(x_norm).to(torch.float32)
        with torch.no_grad():
            y_tensor = model(x_tensor)
        return y_tensor.detach().numpy()

    background = data["x_train"].values
    x_explain = data["x_test"].values
    explainer = shap.KernelExplainer(combined_model, background)
    shap_values = explainer.shap_values(x_explain)
    shap_values = shap_values[0] if isinstance(shap_values, list) else shap_values
    shap_values = np.squeeze(shap_values, axis=2) if shap_values.ndim == 3 else shap_values

    plt.figure()
    shap.summary_plot(shap_values, x_explain, feature_names=x.columns, show=False)
    plt.savefig(result_dir / "shap_summary_plot.png", dpi=300, bbox_inches="tight")
    plt.close()

    df_input = pd.DataFrame(x_explain, columns=x.columns)
    df_shap = pd.DataFrame(shap_values, columns=["shap_" + col for col in x.columns])
    pd.concat([df_input, df_shap], axis=1).to_excel(
        result_dir / "shap_data.xlsx", index=False
    )
    print(f"Saved SHAP outputs for {run_name} to {result_dir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else RUN_ID)
