Uncertainty pipeline
====================

Run the scripts from this directory:
  cd "MachineLearning/Theory-experiment_hybrid_model/NN/uncertainty/code"

Edit parameters in:
  00_config.py

RUN_ID is a neutral local run label. Lightweight public results are written
directly to ../result/. Model weights are written to
../model_load/uq_reference/ and are excluded from the repository.

Running order
-------------

Step 01: 01_train_cv_seed_ensemble.py
  Trains the repeated cross-validation seed ensemble and calculates ML
  epistemic uncertainty (sigma_ML_epi). The ensemble prediction center is not
  exported as a separate public grid field because the central prediction is
  read from MachineLearning/Theory-experiment_hybrid_model/NN/predict/results/analysis.xlsx.

Step 02: 02_build_applicability_domain.py
  Builds the ILR-kNN applicability domain and labels grid points as AD-in or
  AD-out.

Step 03: 03_build_oof_risk_coverage_calibration.py
  Calibrates tau_epi from AD-in out-of-fold risk-coverage. AD status and
  tau_epi divide the grid into Z1-Z4. These zones describe applicability and
  ML epistemic uncertainty; they do not alone qualify a candidate.

Step 04: 04_propagate_descriptor_uncertainty.py
  Passes complete four-descriptor vectors for the evaluated theoretical
  settings through the trained model ensemble. Each source is summarized by
  the standard deviation of its predicted potentials, and the five source
  variances are combined in quadrature to obtain sigma_descriptor.

  Inputs:
    ../input/md_all_compositions_properties.xlsx
    ../input/single_site_results_Co_Ni-Co-Mn-Mn-Mn.xlsx
    ../input/single_site_results_Fe_Co-Zn-Fe-Ni-Mn.xlsx
    ../model_load/uq_reference/cv_seed_ensemble/

  This is an empirical propagation over evaluated settings. It does not use
  Monte Carlo descriptor perturbations or an assumed descriptor distribution.

Step 05: 05_build_uncertainty_table.py
  Estimates the empirical residual term from out-of-fold residuals and repeat
  experiments, then constructs the total uncertainty and formal 95% empirical
  prediction interval:

    sigma_total^2 = sigma_ML_epi^2 + sigma_descriptor^2
                    + sigma_emp_residual^2
    H95 = q0.95 * sigma_total
    PI95_lower = predicted_potential - H95
    PI95_upper = predicted_potential + H95

  Required full-grid output:
    ../result/uncertainty_table.xlsx

Step 06: 06_plot_reliability_maps.py
  Draws the reliability map from the formal uncertainty table.

Step 07: 07_recommend_high_potential_candidates.py
  Produces 40 unique, previously untested experimental candidates through two
  independent channels.

  Global prediction channel (10 points):
    Select the 10 lowest predicted_potential values over the untested grid.
    This channel has no prediction-interval admission gate. Each selected row
    still carries all uncertainty components and the formal PI95.

  Uncertainty-guided channel (30 points):
    Exclude experiments and the global Top-10, require
    PI95_lower < 1.424 V, then select by within-zone optimistic-bound rank and
    ILR composition diversity.

    Z1: 2 zone-top + 4 diversity
    Z2: 3 zone-top + 7 diversity
    Z3: 2 zone-top + 4 diversity
    Z4: 0 zone-top + 8 diversity

  The global Top-10 are included in the ILR distance reference for the 30-point
  channel. If any zone has too few formally eligible points, the script stops;
  candidates that fail the PI criterion are not used to fill the quota.

Step 08: 08_extract_selected_candidate_uncertainty.py
  Reconciles all 40 Step 07 selections against uncertainty_table.xlsx,
  verifies the uncertainty components, interval identities, channel counts,
  and 30-point zone quotas, and writes the reviewer-facing workbook.

Formal result interface
-----------------------

Steps 07-08 require:
  predicted_potential, sigma_ML_epi, sigma_descriptor,
  sigma_emp_residual, sigma_total, H95, PI95_lower, PI95_upper, and Zone.

No temporary interval, estimated half-width, or AD/epistemic-only fallback is
used. If the formal table is missing, run Steps 01-05 first.

Repository outputs
------------------

The repository retains lightweight audit tables, summaries, figures, and the
final 40-point candidate workbooks under ../result/. Full-grid tables and all
model weights remain local and are excluded by .gitignore.

Input and output scope
----------------------

../input/, MachineLearning/data/, and MachineLearning/Theory-experiment_hybrid_model/NN/predict/results/analysis.xlsx are inputs.
The workflow writes only to ../model_load/uq_reference/ and ../result/.
