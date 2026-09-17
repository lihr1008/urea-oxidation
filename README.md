# Urea oxidation workflow

This repository contains the literature-retrieval, evidence-fusion, atomistic-simulation, corpus-analysis, and machine-learning workflows used to study high-entropy Prussian blue analogues for urea oxidation.

## DeepSearcher

`Deepsearcher/` contains the retrieval-augmented literature workflow and its fixed outputs.

- `data/` contains ten Web of Science-derived PBA and UOR corpus partitions used for element screening.
- `query_select_*.py` combines DeepSearcher, BGE-M3 embeddings, and Milvus retrieval with four LLM backends (DeepSeek, GPT, Gemini, and Qwen) to generate candidate element combinations.
- `query_syn_*.py` applies the retrieval workflow to synthesis-oriented literature for protocol recommendation.
- `result/elements/` and `result/synthesis/` contain the DeepSearcher-assisted outputs.
- `without_deepsearcher/` contains the no-retrieval controls used for comparison and stability analysis.

The underlying framework and its dependency specification are included in the same directory. See [Deepsearcher/README.md](Deepsearcher/README.md) for details.

## Knowledge fusion

`KnowledgeFusion/` converts the fixed DeepSearcher element-screening outputs into auditable candidate decisions. The workflow includes input-hash verification, a model-neutral EvidenceAudit interface, anonymized frozen claim and evidence ledgers, reliability discounting, basic probability assignment construction, Yager fusion, a Dempster sensitivity comparison, four non-compensatory evidence gates, and parameter-sensitivity analysis. Numerical results are supplied as CSV/JSON files and a reader-oriented Excel workbook.

See [KnowledgeFusion/README.md](KnowledgeFusion/README.md) and [the EvidenceAudit README](KnowledgeFusion/evidence_audit/README.md) for details.

## NLP corpus analysis

`NLP/` cleans Web of Science records, extracts chemical entities, counts metal occurrences, and constructs element co-occurrence relationships. The supplied PBA corpus includes processed element-frequency and element-link tables. See [NLP/README.md](NLP/README.md) for the processing sequence.

## Atomistic calculations

`AutoCal/DFT/` uses atomate2 and VASP workflows to construct and calculate molecular, slab, and adsorbate-slab systems for different local metal environments. `AutoCal/MD/` generates multimetal initial structures and prepares UFF-based LAMMPS simulations after the GROMACS structure-generation stage.

## MD trajectory preprocessing

`MD preprocess/` contains the trajectory-processing script and batched outputs for the 20,000-composition MD dataset. Production frames are sampled, local six-metal motifs are identified around each metal center, and motif counts are stored for descriptor aggregation.

## MD-derived properties

`MD properties/` combines local-motif statistics with DFT adsorption-energy and Bader-charge data to obtain composition-level theoretical descriptors. It contains the adsorption-energy, CO charge-transfer, and urea charge-transfer processing scripts and their batched results.

## Machine learning

`MachineLearning/` contains the shared data, experiment-only baseline, theory-experiment hybrid models, interpretation analyses, and uncertainty-aware screening.

- `data/` contains the experimental potentials, theory-derived descriptors, and fine-grid composition definitions shared by the workflows.
- `Experiment-only_model/NN/` implements the direct composition-to-`E10` neural-network baseline, including preprocessing, training/test evaluation, cross-validation, fine-grid prediction, element-level SHAP analysis, and ILR-PCA visualization.
- `Theory-experiment_hybrid_model/` implements the two-stage composition-to-four-descriptors-to-`E10` workflow. Its `NN/predict/` directory contains neural-network training, fine-grid screening, element-level SHAP, and hierarchical composition-to-descriptor and descriptor-to-`E10` SHAP analyses. DT, KNN, RF, SVR, and XGBoost implementations provide alternative first- and second-stage regression models.
- `Theory-experiment_hybrid_model/NN/uncertainty/` contains the repeated cross-validation ensemble, ILR-kNN applicability domain, epistemic risk-coverage calibration, five-source descriptor-uncertainty propagation, empirical residual uncertainty, prediction intervals, Z1-Z4 reliability classification, and the 40-composition reliability-aware candidate-selection workflow.

See [MachineLearning/README.md](MachineLearning/README.md) and [the uncertainty pipeline guide](MachineLearning/Theory-experiment_hybrid_model/NN/uncertainty/code/README_pipeline.txt) for the execution order and detailed outputs.
