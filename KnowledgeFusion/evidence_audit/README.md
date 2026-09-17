# EvidenceAudit

This lightweight package documents the semantic evidence-audit interface used before the DS/Yager fusion stage. It prepares deterministic JSONL request batches, validates structured responses returned by an external LLM auditor, and converts accepted audit records into the frozen CSV interface read by `KnowledgeFusion`.

The code is model-client neutral and makes no network or model calls. Request files can be submitted to any LLM service that supports structured JSON output. The returned JSONL must preserve the request and prompt identifiers and hashes.

## Protocol represented here

The 498 frozen claims used by the current analysis were produced by Protocol 7 context-aware atomization. The active sequence is:

1. context-aware claim atomization of each complete DeepSearcher answer;
2. exhaustive claim-to-source mapping within the same DeepSearcher run;
3. controlled evidence classification for mapped claim-source pairs;
4. full same-run rechecking only when a claim is revised during adjudication;
5. human adjudication and conversion to frozen claim/evidence ledgers.

Protocol 7 does not include an independent semantic-merge stage. Claims are not merged across answers or runs.

## Directory contents

- `code/00_prepare_audit_requests.py`: creates atomization, mapping, classification, or revised-claim request batches.
- `code/01_validate_audit_responses.py`: verifies response structure, hashes, IDs, coverage, same-run constraints, controlled vocabularies, and exact excerpt locations.
- `code/02_build_frozen_inputs.py`: applies adjudication decisions and writes `frozen_claims.csv`, `frozen_evidence.csv`, and `frozen_review_summary.json`.
- `code/audit_common.py`: deterministic hashing, request construction, validation, and conversion logic shared by the three entry points.
- `code/test_evidence_audit.py`: end-to-end synthetic test plus read-only validation of the current frozen interface.
- `prompts/`: the four Protocol 7 prompts used by the represented workflow.
- `schemas/`: JSON Schema descriptions for run packets, four response types, and adjudication.
- `examples/`: a small synthetic request/response sequence and its expected frozen outputs. It contains no historical literature excerpts.

## Input contract

`run_packet.schema.json` defines one complete DeepSearcher run. Each JSONL row contains:

- a stable `record_id` and candidate composition;
- the complete frozen final answer used for atomization;
- the source packets recovered from that same run, each with a stable source ID and evidence text;
- optional DOI, UT, publication-year, and declared independent-source metadata. A declared independence group is preserved; DOI/UT/source-packet fallbacks are used only when it is absent.

The same-run source constraint is enforced when mapping and classification requests are generated. A response cannot introduce a claim, source packet, or mapping ID that was absent from its request.

Each generated request records `task`, `request_id`, `request_hash`, `prompt_id`, `prompt_version`, and `prompt_sha256`. The request hash covers the complete canonical request content except the hash field itself. Responses must return the matching request and prompt identifiers and hashes.

## Running the stages

Run commands from `KnowledgeFusion/evidence_audit/code`. Paths below refer to the supplied synthetic example; substitute private run packets and response files when reproducing a full audit.

Prepare atomization requests:

```text
python 00_prepare_audit_requests.py --phase atomization --run-packets ../examples/input/run_packets.jsonl --output atomization_requests.jsonl
```

Submit each request together with the matching file from `../prompts/` to a structured-output LLM auditor. Store one returned JSON object per line, then validate it:

```text
python 01_validate_audit_responses.py --requests atomization_requests.jsonl --responses atomization_responses.jsonl --output validated_atomization.jsonl
```

Use validated atomization output to prepare same-run mapping requests:

```text
python 00_prepare_audit_requests.py --phase mapping --run-packets ../examples/input/run_packets.jsonl --responses validated_atomization.jsonl --output mapping_requests.jsonl
```

After validating the mapping responses, prepare controlled classification requests:

```text
python 00_prepare_audit_requests.py --phase classification --run-packets ../examples/input/run_packets.jsonl --responses validated_mapping.jsonl --output classification_requests.jsonl
```

If adjudication revises a claim, prepare a full same-run recheck:

```text
python 00_prepare_audit_requests.py --phase recheck --run-packets ../examples/input/run_packets.jsonl --revised-claims revised_claims.jsonl --output recheck_requests.jsonl
```

Finally, convert validated records and adjudication decisions to the frozen interface:

```text
python 02_build_frozen_inputs.py --atomization validated_atomization.jsonl --mapping validated_mapping.jsonl --classification validated_classification.jsonl --adjudication adjudication.json --output-dir frozen_output
```

The generated CSV columns match `KnowledgeFusion/input`. `semantic_review_multiplier` is `1` for accepted audit records and the declared semantic discount for machine-unreviewed records. Rejected sources and excluded claims do not enter the DS/Yager inputs.

## Validation

Run the complete synthetic test and the read-only reference check with:

```text
python test_evidence_audit.py
python 02_build_frozen_inputs.py --verify-reference ../../input
```

The validator checks required fields, allowed enumerations, unique IDs, complete request coverage, request and prompt hashes, exact continuous excerpts, one mapping per supplied claim, one revised-claim assessment per same-run source, and classification consistency such as `not_addressed` requiring `entailment=none`.

## Connection to KnowledgeFusion

The output files `frozen_claims.csv`, `frozen_evidence.csv`, and `frozen_review_summary.json` are the semantic-audit boundary. The fusion code reads those files and performs reliability discounting, independent-source collapsing, DS/Yager fusion, non-compensatory gates, and sensitivity analysis. No LLM output is interpreted directly by the fusion code.

## Reproducibility boundary

The public repository retains the frozen DeepSearcher outputs and frozen audit ledgers but not the complete historical source-packet response archive. Therefore, the synthetic example can be run end to end and the preserved 498-claim/2141-evidence interface can be verified, while the historical LLM semantic judgments cannot be regenerated from this repository alone. Missing historical responses are never replaced with approximate data.
