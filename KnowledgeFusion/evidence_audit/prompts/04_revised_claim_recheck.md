Prompt-Version: 7.0.0-revision-recheck
Prompt-ID: revised-claim-full-recheck

# Role

A researcher revised one self-contained claim. Recheck that revised claim against
every supplied source packet from the same original DeepSearcher run. Treat each
source independently.

For each source packet, first return mapping_status as related,
possibly_related, or unrelated with an exact mapping_excerpt and reason. If
related or possibly_related, also return the controlled evidence classification:
relevant, direction, entailment, system_relevance, study_type, traceability,
verified_excerpt, and classification_reason. Use the same definitions as
03_controlled_evidence_classification.md. Do not assign numerical scores.

Return one JSON object with request_id, request_hash, prompt_id, prompt_sha256,
protocol_version, review_model, review_date, and source_assessments. Return
exactly one assessment for each source_packet_id. For unrelated sources, leave
classification fields and excerpts empty.