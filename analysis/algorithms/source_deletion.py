"""Complete fixed posthoc censuses and source-dependency deletions. No tests."""
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import traceback

P = Path(__file__).resolve().parents[1]
B = P/"presentation_core"
S1 = P / "strengthening1_current_contrast_transparency_20261006"
O = P / "strengthening2_all35_landscape_source_influence_20261006"
PROTO = O / "第二项补强_全35基因与全部来源影响_检验前固定事后描述协议.json"

def now():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def json_rows(path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            yield json.loads(line)

def write_json(name, obj):
    (O / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

def write_table(name, rows):
    with (O / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def percent(n, d):
    return 100 * n / d if d else "NA"

def range_record(rows, column, key):
    available = [r for r in rows if isinstance(r[column], (float, int))]
    if not available:
        return {"minimum": "NA", "maximum": "NA", "all_NA": True}
    lo, hi = min(r[column] for r in available), max(r[column] for r in available)
    return {"minimum": lo, "maximum": hi, "minimum_omissions": [r[key] for r in available if r[column] == lo], "maximum_omissions": [r[key] for r in available if r[column] == hi], "available": len(available), "NA": len(rows) - len(available), "meaning": "deterministic full-family range, not CI"}

def main():
    assert not (O / "START_计算.json").exists(), "Inspect old process/receipt before any restart."
    q = json.loads(PROTO.read_text(encoding="utf-8"))
    frozen = json.loads((O / "FREEZE_协议与输入核验.json").read_text(encoding="utf-8"))
    assert sha(PROTO) == frozen["protocol_sha256"]
    for r in q["inputs"]:
        p = Path(r["path"])
        assert p.stat().st_size == r["bytes"] and sha(p) == r["sha256"], r["path"]
    write_json("START_计算.json", {"started_utc": now(), "protocol_sha256": sha(PROTO), "script_sha256": sha(Path(__file__)), "new_tests": 0})
    genes = q["all_candidates"]["genes_alphabetic"]
    sources = q["all_candidates"]["historical_submitter_labels_literal"]
    csources = q["all_candidates"]["current_endpoint_submitter_labels_literal"]
    scopes = q["scope_categories_in_fixed_order"]
    directions = [">".join(d) for d in q["event_directions_in_fixed_order"]]
    gc = {r["gene"]: r for r in read_table(B / "inputs" / "27_all35_gene_descriptive_census.tsv")}
    hc = {r["submitter"]: r for r in read_table(B / "inputs" / "29_all_primary_submitter_concentration.tsv")}
    sc = {r["gene"]: r for r in read_table(B / "inputs" / "33_all35_corrected_scope_counts.tsv")}
    assert sorted(gc) == sorted(sc) == genes and sorted(hc) == sources
    ng, ns, ug, us, eg, es = (Counter() for _ in range(6))
    elig = {}
    seen = set()
    for r in json_rows(B / "inputs" / "40_all_selected_SCV_eligibility_and_source_refs.jsonl.gz"):
        k = r["SCV_base"]
        assert k not in seen
        seen.add(k)
        if r["primary_eligible"]:
            assert r["gene"] in genes and r["submitter"] in sources
            ng[r["gene"]] += 1
            ns[r["submitter"]] += 1
            elig[k] = (r["gene"], r["submitter"], int(r["n_major_updates"]), r["VariationID"])
            if r["n_major_updates"]:
                ug[r["gene"]] += 1
                us[r["submitter"]] += 1
    assert len(seen) == 40834 and len(elig) == 30938
    events = list(json_rows(B / "inputs" / "41_all_primary_major_events.jsonl.gz"))
    assert len(events) == 36
    event_scv = Counter()
    gd, hd, all_d = defaultdict(Counter), defaultdict(Counter), Counter()
    for r in events:
        k, g, s = r["SCV_base"], r["gene"], r["submitter"]
        assert k in elig and elig[k][:2] == (g, s) and str(elig[k][3]) == str(r["VariationID"])
        d = r["left_tier"] + ">" + r["right_tier"]
        assert d in directions and r["right_version"] > r["left_version"]
        event_scv[k] += 1
        eg[g] += 1
        es[s] += 1
        gd[g][d] += 1
        hd[s][d] += 1
        all_d[d] += 1
    assert len(event_scv) == 35
    assert all(v[2] == event_scv[k] for k, v in elig.items())
    for g in genes:
        assert (ng[g], ug[g], eg[g]) == tuple(int(gc[g][x]) for x in ["continuous_primary_SCVs", "SCVs_with_major_update", "major_update_events"]), g
    for s in sources:
        assert (ns[s], us[s], es[s]) == tuple(int(hc[s][x]) for x in ["eligible_SCVs", "updated_distinct_SCVs", "update_events"]), s
    gs, total_scope, scope_ids = defaultdict(Counter), Counter(), set()
    for r in json_rows(B / "inputs" / "47_all4078_original_and_corrected_condition_scope_with_sources.jsonl.gz"):
        vid, g, c = r["VariationID"], r["gene"], r["B_corrected_source_supported_scope"]
        assert vid not in scope_ids and g in genes and c in scopes
        scope_ids.add(vid)
        gs[g][c] += 1
        total_scope[c] += 1
    assert len(scope_ids) == 4078
    for g in genes:
        assert gs[g] == Counter(json.loads(sc[g]["corrected_scope_counts"]))
        assert sum(gs[g].values()) == int(sc[g]["official_PLP_variants"]) == int(gc[g]["latest_official_PLP_variants"])
    current = read_table(S1 / "01_all1141_variant_level_transparency_and_comparability.tsv")
    side_rows = read_table(S1 / "02_all_endpoint_source_records_with_provenance.tsv")
    assert len(current) == 1141 and len(side_rows) == 2320
    vr = {r["VariationID"]: r for r in current}
    assert len(vr) == 1141
    dependencies = defaultdict(set)
    side_count = Counter()
    for r in side_rows:
        vid = r["VariationID"]
        assert vid in vr and r["gene"] == vr[vid]["gene"] and r["Submitter_raw"] in csources
        dependencies[vid].add(r["Submitter_raw"])
        side_count[(vid, r["side"])] += 1
    for vid, r in vr.items():
        assert r["gene"] in genes
        assert side_count[(vid, "earliest")] == int(r["earliest_source_records"])
        assert side_count[(vid, "latest")] == int(r["latest_source_records"])

    def current_stats(records):
        records = list(records)
        return {"contrast_variants": len(records),
                "both_criteria_variants": sum(r["earliest_review_summary"] == r["latest_review_summary"] == "CRITERIA_PROVIDED" for r in records),
                "both_contributing_variants": sum(r["earliest_contribution_summary"] == r["latest_contribution_summary"] == "YES" for r in records),
                "same_condition_variants": sum(r["complete_reported_C_ID_set_comparability"] == "same" for r in records),
                "different_condition_variants": sum(r["complete_reported_C_ID_set_comparability"] == "different" for r in records),
                "unresolved_condition_variants": sum(r["complete_reported_C_ID_set_comparability"] == "unresolved" for r in records),
                "original_window_event_variants": sum(r["has_original_strict_window_update"] == "True" for r in records)}

    base_current = current_stats(current)
    assert list(base_current.values()) == [1141, 900, 893, 232, 659, 250, 19]
    ng_current = {g: current_stats(r for r in current if r["gene"] == g) for g in genes}
    full_gene = []
    full_h = []
    loo_gene = []
    loo_h = []
    loo_c = []
    dir_cells = []
    for g in genes:
        nplp = sum(gs[g].values())
        row = {"gene": g, "stable_endpoint_variants": int(gc[g]["stable_endpoint_variants"]), "eligible_SCVs": ng[g], "updated_distinct_SCVs": ug[g], "major_update_events": eg[g], "observed_updated_SCV_percent": percent(ug[g], ng[g]), "official_PLP_variants": nplp}
        row.update({"scope__" + c: gs[g][c] for c in scopes})
        row.update({"direction__" + d: gd[g][d] for d in directions})
        row.update({"current__" + k: v for k, v in ng_current[g].items()})
        full_gene.append(row)
        kept = {k: base_current[k] - ng_current[g][k] for k in base_current}
        row = {"omitted_gene": g, "removed_eligible_SCVs": ng[g], "removed_updated_SCVs": ug[g], "remaining_eligible_SCVs": 30938 - ng[g], "remaining_updated_SCVs": 35 - ug[g], "remaining_events": 36 - eg[g], "updated_SCV_percent": percent(35 - ug[g], 30938 - ng[g]), "remaining_official_PLP_variants": 4078 - nplp}
        for c in scopes:
            row["remaining_scope__" + c] = total_scope[c] - gs[g][c]
            row["scope_percent__" + c] = percent(total_scope[c] - gs[g][c], 4078 - nplp)
        row.update({"remaining_current__" + k: v for k, v in kept.items()})
        row["current_both_criteria_percent"] = percent(kept["both_criteria_variants"], kept["contrast_variants"])
        row["current_both_contributing_percent"] = percent(kept["both_contributing_variants"], kept["contrast_variants"])
        loo_gene.append(row)
    for s in sources:
        full_h.append({"submitter": s, "eligible_SCVs": ns[s], "updated_distinct_SCVs": us[s], "major_update_events": es[s], "eligible_SCV_share_percent": percent(ns[s], 30938), "updated_SCV_share_percent": percent(us[s], 35), "event_share_percent": percent(es[s], 36)})
        loo_h.append({"omitted_submitter": s, "removed_eligible_SCVs": ns[s], "removed_updated_SCVs": us[s], "removed_events": es[s], "remaining_eligible_SCVs": 30938 - ns[s], "remaining_updated_SCVs": 35 - us[s], "remaining_events": 36 - es[s], "updated_SCV_percent": percent(35 - us[s], 30938 - ns[s])})
    for s in csources:
        removed_ids = {vid for vid in vr if s in dependencies[vid]}
        kept = current_stats(r for r in current if r["VariationID"] not in removed_ids)
        row = {"omitted_endpoint_source_literal": s, "removed_contrast_variants": len(removed_ids), "removed_source_side_rows": sum(r["Submitter_raw"] == s for r in side_rows)}
        row.update({"remaining__" + k: v for k, v in kept.items()})
        row["both_criteria_percent"] = percent(kept["both_criteria_variants"], kept["contrast_variants"])
        row["both_contributing_percent"] = percent(kept["both_contributing_variants"], kept["contrast_variants"])
        row["same_condition_percent"] = percent(kept["same_condition_variants"], kept["contrast_variants"])
        row["different_condition_percent"] = percent(kept["different_condition_variants"], kept["contrast_variants"])
        row["unresolved_condition_percent"] = percent(kept["unresolved_condition_variants"], kept["contrast_variants"])
        loo_c.append(row)
    for level, keys, cmap in [("gene", genes, gd), ("historical_submitter", sources, hd)]:
        for k in keys:
            for d in directions:
                dir_cells.append({"level": level, "omitted_key": k, "direction": d, "source_events": cmap[k][d], "all_events_this_direction": all_d[d], "remaining_events_this_direction": all_d[d] - cmap[k][d]})
    mappings = []
    for prefix, rows, key, coverage in [("D", loo_h, "omitted_submitter", "removed_eligible_SCVs"), ("C", loo_c, "omitted_endpoint_source_literal", "removed_contrast_variants")]:
        for i, r in enumerate(sorted(rows, key=lambda r: (-r[coverage], r[key]))[:2], 1):
            mappings.append({"figure_label": prefix + str(i), "source_literal": r[key], "coverage_metric": coverage, "removed_units": r[coverage], "selection_rule": "largest removed coverage; literal-name tie-break; not sensitivity value or significance"})
    for name, rows in [
        ("01_all35_gene_landscape_corrected_scope_and_separate_units.tsv", full_gene),
        ("02_all265_historical_submitter_census.tsv", full_h),
        ("03_all35_leave_one_gene_descriptors.tsv", loo_gene),
        ("04_all265_leave_one_historical_submitter_descriptors.tsv", loo_h),
        ("05_all96_current_source_dependency_deletions.tsv", loo_c),
        ("06_all1800_gene_source_direction_cells_including_zeros.tsv", dir_cells),
        ("07_fixed_coverage_based_figure_labels.tsv", mappings),
    ]:
        write_table(name, rows)
    summary = {
        "computed_utc": now(), "status": "NUMERIC_COMPLETE_PENDING_INDEPENDENT_CHECK_AND_FIGURE_REVIEW",
        "protocol_sha256": sha(PROTO), "new_tests": 0, "new_models": 0,
        "all_genes": 35, "all_historical_submitter_labels": 265, "all_current_endpoint_source_labels": 96,
        "eligible_SCVs": 30938, "updated_distinct_SCVs": 35, "major_update_events": 36,
        "baseline_observed_updated_SCV_percent": percent(35, 30938),
        "genes_with_updates": sum(ug[g] > 0 for g in genes), "genes_with_zero_updates": sum(ug[g] == 0 for g in genes),
        "sources_with_updates": sum(us[s] > 0 for s in sources), "sources_with_zero_updates": sum(us[s] == 0 for s in sources),
        "all_six_event_directions": {d: all_d[d] for d in directions},
        "official_PLP_variants": 4078, "corrected_scope_counts": {c: total_scope[c] for c in scopes},
        "baseline_current_stats": base_current, "baseline_current_both_criteria_percent": percent(900, 1141),
        "baseline_current_both_contributing_percent": percent(893, 1141),
        "gene_omission_frequency_range": range_record(loo_gene, "updated_SCV_percent", "omitted_gene"),
        "historical_source_omission_frequency_range": range_record(loo_h, "updated_SCV_percent", "omitted_submitter"),
        "gene_omission_current_criteria_range": range_record(loo_gene, "current_both_criteria_percent", "omitted_gene"),
        "current_source_dependency_criteria_range": range_record(loo_c, "both_criteria_percent", "omitted_endpoint_source_literal"),
        "current_source_dependency_contribution_range": range_record(loo_c, "both_contributing_percent", "omitted_endpoint_source_literal"),
        "gene_omission_scope_percent_ranges": {c: range_record(loo_gene, "scope_percent__" + c, "omitted_gene") for c in scopes},
        "coverage_labels": mappings,
        "reconstructed_all_original_gene_submitter_counts_from_existing40_41": True,
        "reconstructed_corrected_all4078_scope_records": True,
        "units_not_patients_and_no_independent_replication": True,
        "original_manuscript_sha256_unchanged": sha(B / "manuscript" / "manuscript_resource_EN_review_ready.md") == frozen["manuscript_sha256"],
    }
    assert summary["original_manuscript_sha256_unchanged"]
    write_json("08_numeric_results_and_bounded_ranges.json", summary)
    files = [{"file": p.name, "bytes": p.stat().st_size, "sha256": sha(p)} for p in sorted(O.glob("0*")) if p.is_file()]
    write_json("COMPLETE_计算输出.json", {"finished_utc": now(), "script_sha256": sha(Path(__file__)), "protocol_sha256": sha(PROTO), "files": files, "clinical_biological_replication": False, "new_tests": 0})
    print(json.dumps({k: v for k, v in summary.items() if k not in ["gene_omission_scope_percent_ranges", "current_source_dependency_contribution_range"]}, ensure_ascii=False))

if __name__ == "__main__":
    try:
        main()
    except Exception:
        if O.exists():
            write_json("FAILURE_计算.json", {"failed_utc": now(), "traceback": traceback.format_exc(), "rerun_authorized_by_this_file": False})
        raise
