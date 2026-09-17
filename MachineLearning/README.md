# Machine-learning workflows

This directory contains source code, shared input data, and selected lightweight results. Trained models, fitted scalers, checkpoints, and run logs are generated locally and are not distributed in the repository.

## Shared data

All workflows resolve their inputs from `MachineLearning/data` by using paths relative to the executing script. Run commands from the corresponding code directory unless a script-specific README states otherwise.

## Experiment-only neural network

The public reference results are stored in `Experiment-only_model/NN/results/experiment_only_reference`.

Run `train.py` first to generate the model and scalers under the ignored local `model_load/experiment_only_reference` directory. Then run `grid_search.py`, `shap_explain.py`, or `dimension_reduction.py` as required. The dimensional-reduction workflow applies ILR before PCA and writes `ilr_pca_results.csv` and `ilr_pca.png` when it is run.

## Theory–experiment hybrid neural network

In `Theory-experiment_hybrid_model/NN/predict/code`, run the scripts in this order:

1. `pretrain.py`
2. `train.py`
3. `grid_search.py`

The first two steps regenerate the model files and scalers required by grid prediction. Lightweight tables and figures remain in `NN/predict/results`; generated model artifacts in that directory are ignored.

The hierarchical SHAP workflow is in `Theory-experiment_hybrid_model/NN/predict/code/shap`. It evaluates composition-to-descriptor and descriptor-to-potential relationships. Run `pretrain.py` and `train.py` first so that the required weights and scalers exist in `NN/predict/results`, then run `shap/run_shap_analysis.py`. Its tables and figures are written to `shap/results`.

The uncertainty workflow is documented separately in `Theory-experiment_hybrid_model/NN/uncertainty/code/README_pipeline.txt`. Its ensemble weights are generated locally and excluded from version control.

## Alternative regression models

For DT, KNN, RF, SVR, and XGBoost, run each `pretrain_*` script before the corresponding `train_*` script. Pretraining creates the first-stage model required by second-stage training. Generated `.pkl` and `.pckl` files are excluded from version control. Where a dimensional-reduction script is provided, it uses the same ILR-PCA reference basis as the neural-network workflow.
