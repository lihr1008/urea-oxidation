from pathlib import Path
import sys

import torch

from nn_common import (
    DirectExperimentModel,
    get_run_dirs_by_id,
    load_pickle,
    read_choose_proportion,
    save_space_prediction,
)

RUN_ID = "experiment_only_reference"

def main(run_id=RUN_ID):
    base_dir = Path(__file__).resolve().parent
    run_name, model_dir, result_dir = get_run_dirs_by_id(base_dir, run_id)
    import json

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
    norm_y = load_pickle(model_dir / "norm_y.pckl")

    x_val = read_choose_proportion()
    x_val_scaled = torch.tensor(norm_x.transform(x_val)).to(torch.float32)
    with torch.no_grad():
        val_pred = model(x_val_scaled).detach().numpy()
    val_pred = norm_y.inverse_transform(val_pred)
    save_space_prediction(result_dir / "space_one_percent.xlsx", x_val, val_pred)
    print(f"Saved space prediction for {run_name} to {result_dir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else RUN_ID)
