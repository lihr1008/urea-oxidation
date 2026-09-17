from pathlib import Path


RUN_ID = "uq_reference"

CODE_DIR = Path(__file__).resolve().parent
BASE_DIR = CODE_DIR.parent
ML_ROOT = Path(__file__).resolve().parents[4]
REPO_ROOT = ML_ROOT.parent

CONFIG = {
    "run_id": RUN_ID,
    "base_seed": 20260326,
    "cv_n_repeats": 5,
    "cv_n_splits": 10,
    "seed_per_fold": 5,
    "prediction_chunk_size": 50000,
    "analysis_path": "MachineLearning/Theory-experiment_hybrid_model/NN/predict/results/analysis.xlsx",
    "repeatability_path": "MachineLearning/Theory-experiment_hybrid_model/NN/uncertainty/input/replicate_experiment.xlsx",
    "repeatability_sheet": "1",
    "md_descriptor_uncertainty_path": "MachineLearning/Theory-experiment_hybrid_model/NN/uncertainty/input/md_all_compositions_properties.xlsx",
    "md_temperature_sheet": "temperature",
    "md_force_field_sheet": "force_field",
    "single_site_descriptor_uncertainty_paths": [
        "MachineLearning/Theory-experiment_hybrid_model/NN/uncertainty/input/single_site_results_Co_Ni-Co-Mn-Mn-Mn.xlsx",
        "MachineLearning/Theory-experiment_hybrid_model/NN/uncertainty/input/single_site_results_Fe_Co-Zn-Fe-Ni-Mn.xlsx",
    ],
    "single_site_functional_sheet": "sheet2_functional",
    "single_site_slab_sheet": "sheet3_big_slab",
    "single_site_second_shell_sheet": "sheet4_second_shell",
    "expected_descriptor_representative_compositions": 5,
    "descriptor_representative_aggregation": "rms",
    "two_point_delta_to_sigma": "sample_sd",
    "descriptor_uncertainty_mode": "global",
    "descriptor_input_space": "physical",
    "candidate_benchmark_E10_V": 1.424,
    "global_top_count": 10,
    "grid_prediction_format": "parquet",
    "ad_quantile": 0.95,
    "ad_k_neighbors": 5,
    "residual_smoothing_k": 5,
    "descriptor_prediction_batch_size": 5000,
    "ilr_zero_replacement_delta": 0.00025,
    "standardized_residual_quantile": 0.95,
    "eps": 1e-12,
    "pretrain": {
        "seed": 0,
        "test_size": 0.2,
        "batch_size": 256,
        "epochs": 500,
        "lr": 0.001,
        "optimizer": "Adam",
        "scheduler_factor": 0.5,
        "scheduler_min_lr": 1e-6,
        "scheduler_patience": 10,
        "dropout": 0.0,
        "weight_decay": 1e-5,
        "input_dim": 5,
        "hidden_units": 512,
        "output_dim": 4,
    },
    "train": {
        "test_size": 0.2,
        "batch_size": 8,
        "epochs": 300,
        "lr": 0.001,
        "optimizer": "Adam",
        "scheduler_factor": 0.3,
        "scheduler_min_lr": 1e-5,
        "scheduler_patience": 10,
        "dropout": 0.4,
        "weight_decay": 1e-5,
        "input_dim": 4,
        "hidden_units": 512,
        "output_dim": 1,
    },
}
