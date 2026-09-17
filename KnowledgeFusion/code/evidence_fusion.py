"""Reliability discounting, Yager fusion, and non-compensatory gating."""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Any, Sequence

from fusion_common import (
    FusionError,
    INPUT_DIR,
    MANIFEST_PATH,
    RESULT_DIR,
    as_bool,
    load_config,
    load_json,
    read_csv,
    text_of,
    write_csv,
)


EXPECTED_SCENARIOS = (
    "nominal",
    "source_strong_discount",
    "source_weak_discount",
    "traceability_strong_discount",
    "traceability_weak_discount",
    "entailment_strong_discount",
    "entailment_weak_discount",
    "gate_strict",
    "gate_lenient",
    "all_strict",
    "all_lenient",
)


def source_strength_key(row: dict[str, Any]) -> str:
    relevance = text_of(row.get("system_relevance"))
    study = text_of(row.get("study_type"))
    if relevance == "direct_PBA_UOR" and study == "experimental":
        return "direct_experimental"
    if relevance == "direct_PBA_UOR" and study in {"computational", "general_chemistry"}:
        return "direct_computational"
    if relevance in {"PBA_only", "UOR_only"} and study == "experimental":
        return "single_domain_experimental"
    if relevance == "substitution":
        return "substitution_adjacent"
    return "foundational_review"


def validate_config(config: dict[str, Any]) -> None:
    profile_specs = (
        ("source_profiles", ("direct_experimental", "direct_computational", "single_domain_experimental", "substitution_adjacent", "foundational_review")),
        ("traceability_profiles", ("exact_chunk_DOI_UT", "abstract_DOI_UT", "identified_record", "partition_only", "untraceable")),
        ("entailment_profiles", ("full", "partial", "contextual", "overreach", "none")),
    )
    for config_key, ordered_categories in profile_specs:
        profiles = config.get(config_key, {})
        if set(profiles) != {"strong_discount", "nominal", "weak_discount"}:
            raise FusionError(f"{config_key} must define strong_discount, nominal, and weak_discount")
        for profile_name, values_by_category in profiles.items():
            if set(values_by_category) != set(ordered_categories):
                raise FusionError(f"{config_key}.{profile_name} has incomplete categories")
            values = [float(values_by_category[key]) for key in ordered_categories]
            if any(value < 0.0 or value > 1.0 for value in values):
                raise FusionError(f"{config_key}.{profile_name} contains a value outside [0,1]")
            if any(left < right for left, right in zip(values, values[1:])):
                raise FusionError(f"{config_key}.{profile_name} must be monotonically non-increasing")

    gates = config.get("gate_profiles", {})
    if set(gates) != {"strict", "nominal", "lenient"}:
        raise FusionError("gate_profiles must define strict, nominal, and lenient")
    required = {"dimension_support", "dimension_oppose", "maximum_ignorance_for_support", "high_conflict"}
    for profile_name, gate in gates.items():
        if set(gate) != required or any(float(gate[key]) < 0.0 or float(gate[key]) > 1.0 for key in required):
            raise FusionError(f"Invalid gate profile: {profile_name}")
    if not (
        gates["strict"]["dimension_support"] >= gates["nominal"]["dimension_support"] >= gates["lenient"]["dimension_support"]
        and gates["strict"]["maximum_ignorance_for_support"] <= gates["nominal"]["maximum_ignorance_for_support"] <= gates["lenient"]["maximum_ignorance_for_support"]
        and gates["strict"]["dimension_oppose"] <= gates["nominal"]["dimension_oppose"] <= gates["lenient"]["dimension_oppose"]
        and gates["strict"]["high_conflict"] <= gates["nominal"]["high_conflict"] <= gates["lenient"]["high_conflict"]
    ):
        raise FusionError("Strict, nominal, and lenient gates are not monotonically ordered")
    scenario_ids = tuple(text_of(row.get("scenario_id")) for row in config.get("sensitivity_scenarios", []))
    if scenario_ids != EXPECTED_SCENARIOS:
        raise FusionError("The 11 preregistered scenarios are missing or out of order")


def scenario_settings(config: dict[str, Any], scenario: dict[str, Any]) -> dict[str, Any]:
    return {
        "scenario_id": scenario["scenario_id"],
        "source_profile": scenario["source_profile"],
        "traceability_profile": scenario["traceability_profile"],
        "entailment_profile": scenario["entailment_profile"],
        "gate_profile": scenario["gate_profile"],
        "source_weights": config["source_profiles"][scenario["source_profile"]],
        "traceability_weights": config["traceability_profiles"][scenario["traceability_profile"]],
        "entailment_weights": config["entailment_profiles"][scenario["entailment_profile"]],
        "gate_thresholds": config["gate_profiles"][scenario["gate_profile"]],
    }


def evidence_alpha(row: dict[str, Any], settings: dict[str, Any]) -> tuple[float, dict[str, float]]:
    source = float(settings["source_weights"].get(source_strength_key(row), 0.0))
    trace = float(settings["traceability_weights"].get(text_of(row.get("traceability")), 0.0))
    entail = float(settings["entailment_weights"].get(text_of(row.get("entailment")), 0.0))
    try:
        semantic = float(row.get("semantic_review_multiplier", 1.0))
    except (TypeError, ValueError):
        semantic = 1.0
    semantic = max(0.0, min(1.0, semantic))
    return source * trace * entail * semantic, {
        "r_source": source,
        "r_traceability": trace,
        "r_entailment": entail,
        "r_semantic": semantic,
    }


def row_bpa(direction: str, alpha: float) -> dict[str, float]:
    alpha = max(0.0, min(1.0, alpha))
    if direction == "support":
        return {"R": alpha, "NR": 0.0, "U": 1.0 - alpha}
    if direction == "oppose":
        return {"R": 0.0, "NR": alpha, "U": 1.0 - alpha}
    if direction == "mixed":
        return {"R": alpha / 2.0, "NR": alpha / 2.0, "U": 1.0 - alpha}
    return {"R": 0.0, "NR": 0.0, "U": 1.0}


def collapse_source(
    rows: Sequence[dict[str, Any]], settings: dict[str, Any]
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    support = oppose = 0.0
    details: list[dict[str, Any]] = []
    for row in rows:
        alpha, factors = evidence_alpha(row, settings)
        direction = text_of(row.get("direction"))
        if direction == "support":
            support = max(support, alpha)
        elif direction == "oppose":
            oppose = max(oppose, alpha)
        elif direction == "mixed":
            support = max(support, alpha / 2.0)
            oppose = max(oppose, alpha / 2.0)
        details.append({**row, **factors, "alpha": alpha})
    if support and oppose:
        effective = max(support, oppose)
        total = support + oppose
        mass = {
            "R": effective * support / total,
            "NR": effective * oppose / total,
            "U": 1.0 - effective,
        }
    elif support:
        mass = row_bpa("support", support)
    elif oppose:
        mass = row_bpa("oppose", oppose)
    else:
        mass = row_bpa("not_addressed", 0.0)
    return mass, details


def _intersect(left: str, right: str) -> str:
    if left == "U":
        return right
    if right == "U":
        return left
    return left if left == right else "EMPTY"


def conjunctive(masses: Sequence[dict[str, float]]) -> dict[str, float]:
    current = {"U": 1.0}
    for mass in masses:
        updated: dict[str, float] = defaultdict(float)
        for left, left_value in current.items():
            for right, right_value in mass.items():
                updated[_intersect(left, right)] += left_value * right_value
        current = dict(updated)
    return {key: current.get(key, 0.0) for key in ("R", "NR", "U", "EMPTY")}


def fuse(masses: Sequence[dict[str, float]], rule: str = "yager") -> dict[str, float]:
    if not masses:
        return {"R": 0.0, "NR": 0.0, "U": 1.0, "K": 0.0, "Bel_R": 0.0, "Pl_R": 1.0, "BetP_R": 0.5}
    raw = conjunctive(masses)
    conflict = raw["EMPTY"]
    if rule == "dempster" and conflict < 1.0 - 1e-12:
        scale = 1.0 / (1.0 - conflict)
        support = raw["R"] * scale
        oppose = raw["NR"] * scale
        uncertainty = raw["U"] * scale
    elif rule == "dempster":
        support, oppose, uncertainty = 0.0, 0.0, 1.0
    else:
        support = raw["R"]
        oppose = raw["NR"]
        uncertainty = raw["U"] + conflict
    total = support + oppose + uncertainty
    support, oppose, uncertainty = support / total, oppose / total, uncertainty / total
    return {
        "R": support,
        "NR": oppose,
        "U": uncertainty,
        "K": conflict,
        "Bel_R": support,
        "Pl_R": support + uncertainty,
        "BetP_R": support + 0.5 * uncertainty,
    }


def dimension_label(row: dict[str, Any], thresholds: dict[str, float]) -> str:
    if row["NR"] >= thresholds["dimension_oppose"]:
        return "opposed"
    if row["K"] >= thresholds["high_conflict"]:
        return "conflict"
    if row["R"] >= thresholds["dimension_support"] and row["U"] <= thresholds["maximum_ignorance_for_support"]:
        return "supported"
    return "insufficient"


def gate_label(dimensions: Sequence[dict[str, Any]]) -> tuple[str, str]:
    labels = [text_of(row.get("dimension_label")) for row in dimensions]
    if "opposed" in labels:
        return "rejected", "at least one necessary design dimension is opposed"
    if "conflict" in labels:
        return "flagged_conflict", "at least one necessary design dimension has high source conflict"
    if labels and all(label == "supported" for label in labels):
        return "supported", "all four necessary design dimensions pass"
    return "flagged_insufficient", "one or more necessary design dimensions remain insufficient"


def evaluate(
    runs: Sequence[dict[str, Any]],
    evidence: Sequence[dict[str, Any]],
    settings: dict[str, Any],
    rule: str,
    critical_dimensions: Sequence[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in evidence:
        if text_of(row.get("included_in_ds")) != "yes":
            continue
        source_id = text_of(row.get("independent_source_id")) or text_of(row.get("source_packet_id"))
        key = (
            text_of(row.get("record_id")),
            text_of(row.get("candidate")),
            text_of(row.get("design_dimension")),
            source_id,
        )
        grouped[key].append(row)

    source_bpas: list[dict[str, Any]] = []
    source_masses: dict[tuple[str, str, str], list[dict[str, float]]] = defaultdict(list)
    for (record_id, candidate, dimension, source_id), rows in grouped.items():
        mass, details = collapse_source(rows, settings)
        source_masses[(record_id, candidate, dimension)].append(mass)
        first = rows[0]
        source_bpas.append({
            "record_id": record_id,
            "candidate": candidate,
            "dimension": dimension,
            "independent_source_id": source_id,
            "source_packet_id": first.get("source_packet_id", ""),
            "DOI": first.get("DOI", ""),
            "UT": first.get("UT", ""),
            "publication_year": first.get("publication_year", ""),
            "scenario_id": settings["scenario_id"],
            "source_profile": settings["source_profile"],
            "traceability_profile": settings["traceability_profile"],
            "entailment_profile": settings["entailment_profile"],
            "gate_profile": settings["gate_profile"],
            "combination_rule": rule,
            "source_m_R": mass["R"],
            "source_m_NR": mass["NR"],
            "source_m_U": mass["U"],
            "collapsed_claim_rows": len(rows),
            "max_alpha": max((item["alpha"] for item in details), default=0.0),
            "minimum_r_semantic": min((item["r_semantic"] for item in details), default=1.0),
            "semantic_review_statuses": "; ".join(
                sorted({text_of(item.get("semantic_review_status")) for item in rows})
            ),
        })

    dimension_rows: list[dict[str, Any]] = []
    decision_rows: list[dict[str, Any]] = []
    for run in runs:
        record_id, candidate = run["record_id"], run["candidate"]
        current: list[dict[str, Any]] = []
        for dimension in critical_dimensions:
            masses = source_masses.get((record_id, candidate, dimension), [])
            result = fuse(masses, rule)
            row = {
                "record_id": record_id,
                "model": run["model"],
                "candidate": candidate,
                "dimension": dimension,
                "scenario_id": settings["scenario_id"],
                "source_profile": settings["source_profile"],
                "traceability_profile": settings["traceability_profile"],
                "entailment_profile": settings["entailment_profile"],
                "gate_profile": settings["gate_profile"],
                "combination_rule": rule,
                **result,
                "independent_source_count": len(masses),
            }
            row["dimension_label"] = dimension_label(row, settings["gate_thresholds"])
            dimension_rows.append(row)
            current.append(row)
        label, reason = gate_label(current)
        decision_rows.append({
            "record_id": record_id,
            "model": run["model"],
            "candidate": candidate,
            "scenario_id": settings["scenario_id"],
            "source_profile": settings["source_profile"],
            "traceability_profile": settings["traceability_profile"],
            "entailment_profile": settings["entailment_profile"],
            "gate_profile": settings["gate_profile"],
            "combination_rule": rule,
            "max_dimension_conflict_K": max((float(row["K"]) for row in current), default=0.0),
            "weakest_dimension_R": min((float(row["R"]) for row in current), default=0.0),
            "maximum_dimension_U": max((float(row["U"]) for row in current), default=1.0),
            "automatic_label": label,
            "gate_reason": reason,
            "supported_dimension_count": sum(row["dimension_label"] == "supported" for row in current),
            "opposed_dimension_count": sum(row["dimension_label"] == "opposed" for row in current),
            "conflict_dimension_count": sum(row["dimension_label"] == "conflict" for row in current),
            "insufficient_dimension_count": sum(row["dimension_label"] == "insufficient" for row in current),
        })
    return dimension_rows, decision_rows, source_bpas


def parameter_candidate_summary(decisions: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in decisions:
        grouped[(text_of(row.get("scenario_id")), text_of(row.get("candidate")))].append(row)
    output: list[dict[str, Any]] = []
    for (scenario_id, candidate), rows in sorted(grouped.items()):
        best = sorted(
            rows,
            key=lambda row: (
                -int(row["supported_dimension_count"]),
                -float(row["weakest_dimension_R"]),
                float(row["max_dimension_conflict_K"]),
                text_of(row.get("record_id")),
            ),
        )[0]
        counts = Counter(text_of(row.get("automatic_label")) for row in rows)
        output.append({
            "scenario_id": scenario_id,
            "candidate": candidate,
            "audited_run_count": len(rows),
            "supported_run_count": counts["supported"],
            "flagged_insufficient_count": counts["flagged_insufficient"],
            "flagged_conflict_count": counts["flagged_conflict"],
            "rejected_count": counts["rejected"],
            "supported_run_rate": counts["supported"] / len(rows) if rows else "",
            "best_record_id": best["record_id"],
            "best_model": best["model"],
            "best_supported_dimension_count": best["supported_dimension_count"],
            "best_weakest_dimension_R": best["weakest_dimension_R"],
            "best_max_dimension_conflict_K": best["max_dimension_conflict_K"],
        })
    return output


def build_run_label_stability(decisions: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in decisions:
        grouped[text_of(row.get("record_id"))].append(row)
    output: list[dict[str, Any]] = []
    for record_id, rows in sorted(grouped.items()):
        by_scenario = {text_of(row.get("scenario_id")): row for row in rows}
        nominal = by_scenario["nominal"]
        changed = [
            scenario_id
            for scenario_id in EXPECTED_SCENARIOS
            if by_scenario[scenario_id]["automatic_label"] != nominal["automatic_label"]
        ]
        drivers: list[str] = []
        if any(item.startswith("source_") for item in changed):
            drivers.append("source_weight")
        if any(item.startswith("traceability_") for item in changed):
            drivers.append("traceability_weight")
        if any(item.startswith("entailment_") for item in changed):
            drivers.append("entailment_weight")
        if any(item.startswith("gate_") for item in changed):
            drivers.append("gate_threshold")
        if any(item.startswith("all_") for item in changed):
            drivers.append("combined_parameters")
        labels = [text_of(by_scenario[item].get("automatic_label")) for item in EXPECTED_SCENARIOS]
        output.append({
            "record_id": record_id,
            "model": nominal["model"],
            "candidate": nominal["candidate"],
            "nominal_label": nominal["automatic_label"],
            "label_set": "; ".join(sorted(set(labels))),
            "scenario_count": len(EXPECTED_SCENARIOS),
            "nominal_label_agreement_count": sum(label == nominal["automatic_label"] for label in labels),
            "nominal_label_agreement_rate": sum(label == nominal["automatic_label"] for label in labels) / len(labels),
            "changed_scenarios": "; ".join(changed),
            "variation_drivers": "; ".join(drivers),
            "parameter_stable": "yes" if not changed else "no",
        })
    return output


def build_rule_sensitivity(
    yager_dimensions: Sequence[dict[str, Any]],
    dempster_dimensions: Sequence[dict[str, Any]],
    yager_decisions: Sequence[dict[str, Any]],
    dempster_decisions: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    dempster_index = {(row["record_id"], row["dimension"]): row for row in dempster_dimensions}
    yager_run = {row["record_id"]: row for row in yager_decisions}
    dempster_run = {row["record_id"]: row for row in dempster_decisions}
    output: list[dict[str, Any]] = []
    for row in yager_dimensions:
        other = dempster_index[(row["record_id"], row["dimension"])]
        output.append({
            "record_id": row["record_id"],
            "model": row["model"],
            "candidate": row["candidate"],
            "dimension": row["dimension"],
            "yager_R": row["R"],
            "yager_NR": row["NR"],
            "yager_U": row["U"],
            "yager_K": row["K"],
            "yager_dimension_label": row["dimension_label"],
            "dempster_R": other["R"],
            "dempster_NR": other["NR"],
            "dempster_U": other["U"],
            "dempster_K": other["K"],
            "dempster_dimension_label": other["dimension_label"],
            "dimension_label_changed": "yes" if row["dimension_label"] != other["dimension_label"] else "no",
            "yager_run_label": yager_run[row["record_id"]]["automatic_label"],
            "dempster_run_label": dempster_run[row["record_id"]]["automatic_label"],
            "run_label_changed": "yes" if yager_run[row["record_id"]]["automatic_label"] != dempster_run[row["record_id"]]["automatic_label"] else "no",
        })
    return output


def _parameter_stability(rows: Sequence[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[text_of(row.get("candidate"))].append(row)
    output: dict[str, dict[str, Any]] = {}
    for candidate, items in grouped.items():
        by_scenario = {text_of(row.get("scenario_id")): row for row in items}
        supported = [int(row["supported_run_count"]) for row in items]
        output[candidate] = {
            "nominal_supported_run_count": int(by_scenario["nominal"]["supported_run_count"]),
            "supported_run_count_min": min(supported),
            "supported_run_count_max": max(supported),
            "scenarios_with_any_supported_run": sum(value > 0 for value in supported),
            "parameter_scenario_count": len(items),
            "supported_count_stable_across_parameters": "yes" if len(set(supported)) == 1 else "no",
        }
    return output


def build_candidate_summary(
    manifest: Sequence[dict[str, Any]],
    decisions: Sequence[dict[str, Any]],
    dimensions: Sequence[dict[str, Any]],
    parameter_candidates: Sequence[dict[str, Any]],
    evidence: Sequence[dict[str, Any]],
    critical_dimensions: Sequence[str],
    semantic_discount: float,
) -> list[dict[str, Any]]:
    assisted = [row for row in manifest if text_of(row.get("source_group")) == "with_deepsearcher"]
    controls = [row for row in manifest if text_of(row.get("source_group")) == "without_deepsearcher"]
    generation_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    control_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    decision_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    dimension_index = {(row["record_id"], row["dimension"]): row for row in dimensions}
    for row in assisted:
        generation_groups[text_of(row.get("candidate"))].append(row)
    for row in controls:
        control_groups[text_of(row.get("candidate"))].append(row)
    for row in decisions:
        decision_groups[text_of(row.get("candidate"))].append(row)

    source_groups: dict[str, dict[str, set[tuple[str, str]]]] = defaultdict(
        lambda: {"reviewed": set(), "auditor_only": set()}
    )
    overreach_counts: Counter[str] = Counter()
    for row in evidence:
        candidate = text_of(row.get("candidate"))
        if text_of(row.get("entailment")) == "overreach":
            overreach_counts[candidate] += 1
        if text_of(row.get("included_in_ds")) != "yes":
            continue
        key = (text_of(row.get("record_id")), text_of(row.get("source_packet_id")))
        if text_of(row.get("semantic_review_status")) in {"human_accepted", "human_revised"}:
            source_groups[candidate]["reviewed"].add(key)
        elif text_of(row.get("semantic_review_status")) == "machine_unreviewed":
            source_groups[candidate]["auditor_only"].add(key)

    stability = _parameter_stability(parameter_candidates)
    output: list[dict[str, Any]] = []
    for candidate, rows in sorted(decision_groups.items()):
        best = sorted(
            rows,
            key=lambda row: (
                -int(row["supported_dimension_count"]),
                -float(row["weakest_dimension_R"]),
                float(row["max_dimension_conflict_K"]),
                text_of(row.get("record_id")),
            ),
        )[0]
        counts = Counter(text_of(row.get("automatic_label")) for row in rows)
        generated = generation_groups[candidate]
        control_count = len(control_groups.get(candidate, []))
        result = {
            "candidate": candidate,
            "generation_count": len(generated),
            "control_generation_count": control_count,
            "generation_count_difference": len(generated) - control_count,
            "models_with_generation": "; ".join(sorted({text_of(row.get("model")) for row in generated})),
            "audited_run_count": len(rows),
            "supported_run_count": counts["supported"],
            "flagged_insufficient_count": counts["flagged_insufficient"],
            "flagged_conflict_count": counts["flagged_conflict"],
            "rejected_count": counts["rejected"],
            "best_record_id": best["record_id"],
            "best_model": best["model"],
            "best_automatic_label": best["automatic_label"],
            "best_supported_dimension_count": best["supported_dimension_count"],
            "best_weakest_dimension_R": best["weakest_dimension_R"],
            "best_max_dimension_conflict_K": best["max_dimension_conflict_K"],
            **stability[candidate],
            "reviewed_source_count": len(source_groups[candidate]["reviewed"]),
            "auditor_only_source_count": len(source_groups[candidate]["auditor_only"]),
            "auditor_overreach_count": overreach_counts[candidate],
            "semantic_discount_for_unreviewed_sources": semantic_discount,
            "evidence_based_assessment": "priority candidate" if counts["supported"] else "insufficient-evidence alternative",
        }
        for dimension in critical_dimensions:
            item = dimension_index[(best["record_id"], dimension)]
            result[f"{dimension}_R"] = item["R"]
            result[f"{dimension}_NR"] = item["NR"]
            result[f"{dimension}_U"] = item["U"]
            result[f"{dimension}_K"] = item["K"]
            result[f"{dimension}_label"] = item["dimension_label"]
        output.append(result)
    return sorted(
        output,
        key=lambda row: (-int(row["supported_run_count"]), -int(row["generation_count"]), row["candidate"]),
    )


def run_pipeline(write_results: bool = True) -> dict[str, list[dict[str, Any]]]:
    config = load_config()
    validate_config(config)
    manifest = load_json(MANIFEST_PATH)
    review = load_json(INPUT_DIR / "frozen_review_summary.json")
    claims = read_csv(INPUT_DIR / "frozen_claims.csv")
    evidence = read_csv(INPUT_DIR / "frozen_evidence.csv")
    if not as_bool(review.get("ready_for_fusion")) or int(review.get("pending_actions", -1)) != 0:
        raise FusionError("The frozen semantic review is not complete")

    runs = [
        {
            "record_id": text_of(row.get("record_id")),
            "model": text_of(row.get("model")),
            "candidate": text_of(row.get("candidate")),
        }
        for row in manifest.get("runs", [])
        if text_of(row.get("source_group")) == "with_deepsearcher"
    ]
    if len(runs) != 40:
        raise FusionError(f"Knowledge fusion requires 40 assisted runs; found {len(runs)}")
    run_ids = {row["record_id"] for row in runs}
    if any(text_of(row.get("record_id")) not in run_ids for row in claims + evidence):
        raise FusionError("A frozen claim or evidence row refers to an unknown assisted run")

    critical_dimensions = config["critical_dimensions"]
    scenario_results: dict[str, tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]] = {}
    parameter_decisions: list[dict[str, Any]] = []
    for scenario in config["sensitivity_scenarios"]:
        settings = scenario_settings(config, scenario)
        dimensions, decisions, source_bpas = evaluate(
            runs, evidence, settings, config["primary_combination_rule"], critical_dimensions
        )
        scenario_results[settings["scenario_id"]] = (dimensions, decisions, source_bpas)
        parameter_decisions.extend(dict(row) for row in decisions)

    nominal_dimensions, nominal_decisions, nominal_bpas = scenario_results["nominal"]
    nominal_settings = scenario_settings(config, config["sensitivity_scenarios"][0])
    dempster_dimensions, dempster_decisions, _ = evaluate(
        runs, evidence, nominal_settings, config["secondary_combination_rule"], critical_dimensions
    )
    human_evidence = [
        row
        for row in evidence
        if text_of(row.get("included_in_ds")) == "yes"
        and text_of(row.get("semantic_review_status")) in {"human_accepted", "human_revised"}
    ]
    human_dimensions, human_decisions, _ = evaluate(
        runs, human_evidence, nominal_settings, config["primary_combination_rule"], critical_dimensions
    )

    strict_supported = {
        row["record_id"]
        for row in scenario_results["gate_strict"][1]
        if row["automatic_label"] == "supported"
    }
    nominal_supported = {
        row["record_id"] for row in nominal_decisions if row["automatic_label"] == "supported"
    }
    lenient_supported = {
        row["record_id"]
        for row in scenario_results["gate_lenient"][1]
        if row["automatic_label"] == "supported"
    }
    if not strict_supported <= nominal_supported <= lenient_supported:
        raise FusionError("Gate-only supported sets violate strict subset nominal subset lenient")

    human_dimension_index = {(row["record_id"], row["dimension"]): row for row in human_dimensions}
    human_decision_index = {row["record_id"]: row for row in human_decisions}
    dependency_by_run: Counter[str] = Counter()
    human_comparison: list[dict[str, Any]] = []
    for row in nominal_dimensions:
        human = human_dimension_index[(row["record_id"], row["dimension"])]
        dependent = row["dimension_label"] != human["dimension_label"]
        if dependent:
            dependency_by_run[row["record_id"]] += 1
        human_comparison.append({
            "record_id": row["record_id"],
            "model": row["model"],
            "candidate": row["candidate"],
            "dimension": row["dimension"],
            "full_evidence_label": row["dimension_label"],
            "human_only_label": human["dimension_label"],
            "auditor_dependency_flag": "yes" if dependent else "no",
            "full_evidence_R": row["R"],
            "full_evidence_NR": row["NR"],
            "full_evidence_U": row["U"],
            "full_evidence_K": row["K"],
            "human_only_R": human["R"],
            "human_only_NR": human["NR"],
            "human_only_U": human["U"],
            "human_only_K": human["K"],
        })
    for row in nominal_decisions:
        human = human_decision_index[row["record_id"]]
        dependent_count = dependency_by_run[row["record_id"]]
        row["human_only_automatic_label"] = human["automatic_label"]
        row["auditor_dependent_dimension_count"] = dependent_count
        row["auditor_dependency_flag"] = "yes" if dependent_count else "no"

    parameter_candidates = parameter_candidate_summary(parameter_decisions)
    run_stability = build_run_label_stability(parameter_decisions)
    stability_index = {row["record_id"]: row for row in run_stability}
    for row in nominal_decisions:
        row["parameter_stable"] = stability_index[row["record_id"]]["parameter_stable"]
        row["nominal_label_agreement_rate"] = stability_index[row["record_id"]]["nominal_label_agreement_rate"]

    semantic_discount = float(review["semantic_discount_for_unreviewed_sources"])
    candidate_summary = build_candidate_summary(
        manifest.get("runs", []),
        nominal_decisions,
        nominal_dimensions,
        parameter_candidates,
        evidence,
        critical_dimensions,
        semantic_discount,
    )
    rule_sensitivity = build_rule_sensitivity(
        nominal_dimensions, dempster_dimensions, nominal_decisions, dempster_decisions
    )

    for row in nominal_dimensions:
        if abs(float(row["R"]) + float(row["NR"]) + float(row["U"]) - 1.0) > 1e-12:
            raise FusionError(f"Dimension masses do not sum to one: {row['record_id']} {row['dimension']}")
    for row in nominal_bpas:
        if abs(float(row["source_m_R"]) + float(row["source_m_NR"]) + float(row["source_m_U"]) - 1.0) > 1e-12:
            raise FusionError(f"Source BPA does not sum to one: {row['record_id']} {row['independent_source_id']}")

    expected = config["reference_expectations"]
    if len(nominal_decisions) != int(expected["assisted_run_count"]):
        raise FusionError("Nominal run-decision count does not match the frozen reference")
    if len(nominal_dimensions) != int(expected["nominal_dimension_count"]):
        raise FusionError("Nominal dimension count does not match the frozen reference")
    if len(config["sensitivity_scenarios"]) != int(expected["scenario_count"]):
        raise FusionError("Sensitivity-scenario count does not match the frozen reference")
    if len(candidate_summary) != int(expected["candidate_count"]):
        raise FusionError("Candidate count does not match the frozen reference")
    priority = next(
        (row for row in candidate_summary if row["candidate"] == expected["priority_candidate"]),
        None,
    )
    if priority is None:
        raise FusionError("The frozen priority candidate is absent")
    if (
        int(priority["generation_count"]) != int(expected["priority_generation_count"])
        or int(priority["supported_run_count"]) != int(expected["priority_supported_run_count"])
        or priority["evidence_based_assessment"] != "priority candidate"
    ):
        raise FusionError("The priority-candidate reference checks failed")

    outputs = {
        "candidate_summary": candidate_summary,
        "run_candidate_decisions": nominal_decisions,
        "run_dimension_masses": nominal_dimensions,
        "source_bpas": nominal_bpas,
        "parameter_sensitivity": parameter_candidates,
        "run_label_stability": run_stability,
        "rule_sensitivity": rule_sensitivity,
        "human_only_comparison": human_comparison,
    }
    if write_results:
        for name, rows in outputs.items():
            write_csv(RESULT_DIR / f"{name}.csv", rows)
    return outputs
