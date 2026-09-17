# KnowledgeFusion

This directory contains the lightweight, reproducible evidence-audit and knowledge-fusion stages used to prioritize multimetal PBA/UOR candidates. The public sequence is input verification → neutral evidence-audit interface → frozen claim/evidence inputs → DS/Yager fusion. The workflow reads the fixed DeepSearcher outputs already stored in this repository and does not modify them.

## Inputs

- `Deepsearcher/result/elements`: 40 DeepSearcher-assisted element-screening outputs.
- `Deepsearcher/without_deepsearcher`: 40 no-retrieval control outputs used only for generation-stability statistics.
- `Deepsearcher/data`: 10 PBA/UOR corpus partitions used for source-consistency verification.
- `KnowledgeFusion/input`: anonymized frozen claim and evidence ledgers required by the fusion calculation.
- `KnowledgeFusion/evidence_audit`: Protocol 7 prompts, request/response schemas, neutral batch utilities, and a synthetic end-to-end example for producing the frozen ledgers.

`input/deepsearcher_manifest.json` binds every fixed source file and frozen ledger to a SHA-256 digest. Verification fails on a missing, additional, renamed, or modified input; it does not fall back to historical paths or approximate matches.

## Method

The current 498 claims come from Protocol 7 context-aware atomization, followed by same-run claim-source mapping, controlled evidence classification, and human adjudication. This workflow has no separate semantic-merge stage. The public `evidence_audit` package does not call an LLM; it prepares auditable requests, validates externally supplied structured responses, and converts accepted records to the frozen input interface.

Each included evidence item is reliability-discounted as

`alpha = source relevance × traceability × entailment × semantic review multiplier`.

Duplicate claims from the same independent source are collapsed before fusion. Yager's rule is primary; Dempster's rule is only a sensitivity comparison. Four necessary dimensions—UOR activity, PBA structural compatibility, wet-chemical synthesis, and alkaline stability/corrosion—are gated independently and cannot compensate for one another. Eleven preregistered scenarios test discount and gate sensitivity.

The masses `R`, `NR`, and `U` express evidential support, opposition, and uncertainty for experimental prioritization; they are not probabilities of experimental success.

## Run

From `KnowledgeFusion/code`:

```text
python 00_verify_deepsearcher_inputs.py
python 01_run_knowledge_fusion.py
python test_knowledge_fusion.py
```

The first command verifies the fixed inputs. The second repeats verification and regenerates the CSV/JSON numerical results. The supplied Excel workbook is a reader-oriented mirror of those machine-readable results and contains no additional calculation logic.

The evidence-audit interface and its synthetic example can be tested separately from `KnowledgeFusion/evidence_audit/code`:

```text
python test_evidence_audit.py
python 02_build_frozen_inputs.py --verify-reference ../../input
```

See `evidence_audit/README.md` for the staged request, validation, adjudication, and frozen-input commands.

## Output

- `candidate_summary.csv`: candidate-level generation, gate, and stability summary.
- `run_candidate_decisions.csv`: 40 nominal run decisions.
- `run_dimension_masses.csv`: 160 nominal dimension-level masses.
- `source_bpas.csv`: nominal independent-source basic probability assignments.
- `parameter_sensitivity.csv`: candidate results for all 11 parameter scenarios.
- `run_label_stability.csv`: run-level sensitivity across those scenarios.
- `rule_sensitivity.csv`: nominal Yager versus Dempster comparison.
- `human_only_comparison.csv`: full-evidence versus reviewed-only comparison.
- `knowledge_fusion_results.xlsx`: the same results arranged for reading.

Under the nominal setting, `Fe-Co-Ni-Mn-Zn` remains the priority candidate: it appears in 27 of the 40 assisted runs, and 3 runs pass all four necessary-dimension gates.

## Scope limitation

The repository does not contain the full historical retrieval-fragment logs, local retrieval database, or complete historical LLM source-packet response archive. The preserved DeepSearcher outputs and frozen audit inputs can be verified, and the DS/Yager fusion can be fully recomputed. The supplied synthetic audit example is also complete, but the historical semantic judgments and earlier evidence-fragment recovery stage cannot be recreated from zero. No substitute data are used when those historical traces are absent.
