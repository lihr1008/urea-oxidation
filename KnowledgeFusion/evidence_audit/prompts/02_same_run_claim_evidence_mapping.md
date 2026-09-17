Prompt-Version: 7.0.0-mapping
Prompt-ID: exhaustive-same-run-mapping

# Role

You are performing a high-recall topic mapping between self-contained claims and
one source segment recovered from the same original DeepSearcher run.

# Task

Read the entire evidence_text once, then assess every supplied claim. Return only
whether the source addresses the scientific subject of the claim:

- related: the source directly discusses the same subject, object, property, or
  mechanism.
- possibly_related: the source discusses a close but incomplete, adjacent, or
  more general version of the subject.
- unrelated: it does not address the claim.

This phase is only a map. Do not decide support, opposition, entailment,
overreach, evidence quality, reliability, candidate value, or probability.
Do not use model identity, recommendation frequency, NLP frequency, retrieval
score, or expected conclusions.

For related and possibly_related, copy a short continuous mapping_excerpt exactly
from evidence_text. For unrelated, mapping_excerpt must be empty. Always give a
brief mapping_reason.

# Output

Return one JSON object with request_id, request_hash, prompt_id, prompt_sha256,
protocol_version, review_model, review_date, and mappings. Each mapping item
contains claim_id, mapping_status, mapping_excerpt, and mapping_reason.
Return exactly one item for every supplied claim_id.
