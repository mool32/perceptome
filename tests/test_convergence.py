"""Convergence subpackage — bundled paper reference, assignment, within-patient test."""

import numpy as np
import pandas as pd
import pytest

import perceptome as pct
from perceptome import convergence as cv

REF = cv.load_convergence_reference()
S1, S2, S4, S6, S7, S8 = (REF["tables"][k] for k in ("S1", "S2", "S4", "S6", "S7", "S8"))


# ── bundled Supporting-Information tables ────────────────────────────────────

def test_tables_have_the_paper_dimensions():
    assert [len(t) for t in (S1, S2, S4, S6, S7, S8)] == [25, 154, 44, 154, 44, 25]


def test_table_names_match_the_bundled_reference():
    hpa = pct.load_hpa_perceptivity()
    assert set(S2.index) == set(S6.index) == set(hpa["R"].index)
    assert set(S4.index) == set(S7.index) == set(pct.list_modules())
    assert set(S8["cell_of_origin"]) <= set(hpa["R"].index)


def test_paper_stated_totals():
    assert S1["n_cells_QC"].sum() == 448399                  # 'full quality-controlled pool ... 448,399 cells'
    assert S6["n_cells_assigned"].sum() == 574451            # 'Table S6 tabulates 574,451 assignments'
    assert S4["mean_contribution"].sum() == pytest.approx(0.711, abs=0.002)   # 'the 44 contributions sum to ... 0.711'
    assert (S7["core_suppression"] + S7["target_distinctiveness"]).abs().max() < 1e-9


def test_table_s2_is_a_distribution_over_152_candidates():
    cand = S2[S2["qc_status"] == "candidate"]
    assert len(cand) == 152
    assert cand["occupancy_pct"].sum() == pytest.approx(100.0, abs=0.5)
    assert set(S2.index[S2["qc_status"] == "quarantined"]) == {"erythrocytes", "platelets"}
    assert cv.quarantined_cell_types() == ["erythrocytes", "platelets"]


def test_cytotrophoblast_and_megakaryocyte_lead_recurrence_and_depth():
    cand = S2[S2["qc_status"] == "candidate"].sort_values("recurrence", ascending=False)
    assert set(cand.index[:2]) == {"cytotrophoblasts", "megakaryocytes"}
    assert (cand["recurrence"].iloc[:2] == 18).all()
    assert cand["recurrence"].iloc[2] == 13 and cand.index[2] == "gastric progenitor cells"
    # and the only pair leading median depth
    assert set(cand.sort_values("median_depth_pct", ascending=False).index[:2]) == {"cytotrophoblasts", "megakaryocytes"}


def test_region_occupancies_reproduce_the_published_enrichment():
    """Fig. 3F: 23.9% / 41.9% / 76.9% of malignant cells; 18.2x / 15.9x / 8.4x over chance."""
    occ = S2["occupancy_pct"]
    pub = {e["boundary"]: e for e in REF["headline"]["enrichment"]}
    for name in ("anchor", "core", "neighbourhood"):
        members = cv.region_members(name)
        pct_in = occ[members].sum()
        assert pct_in == pytest.approx(pub[name]["occupancy_pct"], abs=0.05), name
        assert cv.enrichment(pct_in / 100, len(members)) == pytest.approx(pub[name]["fold"], abs=0.1), name
    members = cv.region_members("region")
    assert len(members) == 14
    assert occ[members].sum() == pytest.approx(pub["region"]["occupancy_pct"], abs=1.0)   # places 13-14 are a tie
    assert cv.enrichment(occ[members].sum() / 100, 14) == pytest.approx(pub["region"]["fold"], abs=0.15)


def test_published_enrichment_arithmetic():
    for e in REF["headline"]["enrichment"]:
        assert e["chance_pct"] == pytest.approx(100 * e["k"] / 152, abs=0.06)
        assert e["occupancy_pct"] / e["chance_pct"] == pytest.approx(e["fold"], rel=0.03)


def test_core_cells_pass_and_low_complexity_targets_are_flagged_in_table_s6():
    """Note S3: core-assigned cells carry 1,396 (megakaryocyte) and 2,430 (cytotrophoblast) median genes;
    the quarantined erythrocyte sits at 890, below the 1,000-gene floor; Note S-Q names three more."""
    floor = REF["definitions"]["qc"]["min_genes_per_cell"]
    med = S6["median_ngenes"]
    assert med["megakaryocytes"] == 1396 and med["cytotrophoblasts"] == 2430
    assert med["erythrocytes"] == 890 < floor
    assert (med[["late spermatids", "gastric chief cells", "pancreatic acinar cells"]] < floor).all()
    assert (med[["cytotrophoblasts", "megakaryocytes"]] > floor).all()


def test_within_patient_registry_is_internally_consistent():
    pw = cv.published_within_patient()
    four = ["colorectal_pelka", "colorectal_hta", "kidney_ccrcc", "liver_hcc_sun"]
    pooled = pw.loc["pooled_four_datasets"]
    assert pw.loc[four, "n_deeper"].sum() == pooled["n_deeper"] == 77
    assert pw.loc[four, "n_patients"].sum() == pooled["n_patients"] == 86
    counted = pw.dropna(subset=["n_deeper", "n_patients"])
    assert (counted["n_deeper"] <= counted["n_patients"]).all()
    assert pw.loc["colorectal_pelka", ["normal_pct", "malignant_pct"]].tolist() == [3.2, 18.9]
    assert pw.loc["lung_gse131907", "n_deeper"] == pw.loc["lung_gse131907", "n_patients"] == 10


def test_registry_carries_the_retractions_and_limits():
    retracted = " ".join(REF["retracted"])
    for phrase in ("Baseline gating", "fold-change", "Stayton", "polyploidy", "cGAS-STING"):
        assert phrase in retracted
    limits = " ".join(REF["limits"])
    assert "not a proven dynamical attractor" in limits and "catalog" in limits


def test_region_members_validation():
    assert cv.region_members("core") == ["cytotrophoblasts", "megakaryocytes"]
    assert set(cv.region_members("core")) < set(cv.region_members("neighbourhood"))
    assert cv.region_members(3) == ["cytotrophoblasts", "megakaryocytes", "gastric progenitor cells"]
    assert cv.region_members(["hepatocytes"]) == ["hepatocytes"]
    with pytest.raises(ValueError):
        cv.region_members("bogus")
    with pytest.raises(ValueError):
        cv.region_members(["not a cell type"])


# ── assignment ───────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def R():
    return pct.load_hpa_perceptivity()["R"].astype(float)


def test_reference_cell_types_assign_to_themselves(R):
    cand = R.drop(index=cv.quarantined_cell_types())
    out = cv.assign_to_reference(cand, warn_scale=False)
    assert (out["nearest"] == cand.index).all()
    assert (out["cosine"] > 0.999).all()
    assert not out["displaced"].any()
    assert out.attrs["calibration"] == "pseudobulk"


def test_quarantined_cell_types_are_reassigned_to_a_surviving_neighbour(R):
    q = R.loc[cv.quarantined_cell_types()]
    out = cv.assign_to_reference(q, warn_scale=False)
    assert (out["first_choice"] == q.index).all() and out["displaced"].all()
    assert not out["nearest"].isin(cv.quarantined_cell_types()).any()
    kept = cv.assign_to_reference(q, quarantine=False, warn_scale=False)
    assert (kept["nearest"] == q.index).all()


def test_shape_only_assignment_ignores_overall_magnitude(R):
    cand = R.drop(index=cv.quarantined_cell_types())
    shifted = cand + 3.0
    shape = cv.assign_to_reference(cand, shape_only=True, warn_scale=False)["nearest"]
    shape_shifted = cv.assign_to_reference(shifted, shape_only=True, warn_scale=False)["nearest"]
    assert (shape == shape_shifted).all()
    raw = cv.assign_to_reference(cand, warn_scale=False)["nearest"]
    raw_shifted = cv.assign_to_reference(shifted, warn_scale=False)["nearest"]
    assert (raw != raw_shifted).mean() > 0.1      # magnitude does move raw assignments


def test_per_cell_scale_scores_warn_in_assign_to_reference(R):
    with pytest.warns(UserWarning, match="single_cell_scaled"):
        cv.assign_to_reference(R.iloc[:30] / 10.0)


def test_the_single_cell_calibration_is_bundled_and_its_absence_points_to_the_build_script(monkeypatch, tmp_path):
    assert pct.load_calibration("single_cell_scaled").name == "single_cell_scaled"
    from perceptome.eigenspace import calibration as calmod
    monkeypatch.setattr(calmod, "_SINGLE_CELL_FILE", tmp_path / "missing.json")
    with pytest.raises(FileNotFoundError, match="13_build_single_cell_calibration"):
        pct.load_calibration("single_cell_scaled")


def test_core_fraction_region_occupancy_and_enrichment():
    near = pd.Series(["cytotrophoblasts"] * 2 + ["megakaryocytes"] * 1 + ["hepatocytes"] * 7)
    assert cv.core_fraction(near) == pytest.approx(0.3)
    assert cv.core_fraction(near, "anchor") == pytest.approx(0.2)
    occ = cv.region_occupancy(near, regions=("anchor", "core"))
    assert occ.loc["core", "occupancy_pct"] == pytest.approx(30.0)
    assert occ.loc["core", "fold_over_chance"] == pytest.approx(0.3 / (2 / 152))
    assert cv.enrichment(0.239, 2) == pytest.approx(18.2, abs=0.05)
    assert cv.enrichment(0.111, 1) == pytest.approx(16.9, abs=0.05)


# ── within-patient test ──────────────────────────────────────────────────────

def _cohort(n_patients, malignant_hits, normal_hits, n_cells=40):
    """Synthetic paired cohort: per patient, `n_cells` cells per arm, given core hits per arm."""
    labels, patient, arm = [], [], []
    for p in range(n_patients):
        for a, hits in (("tumour", malignant_hits[p]), ("normal", normal_hits[p])):
            labels += ["cytotrophoblasts"] * hits + ["hepatocytes"] * (n_cells - hits)
            patient += [f"P{p:02d}"] * n_cells
            arm += [a] * n_cells
    idx = pd.RangeIndex(len(labels))
    return pd.Series(labels, index=idx), pd.Series(patient, index=idx), pd.Series(arm, index=idx)


def test_within_patient_reproduces_the_papers_exact_statistics():
    """10 of 10 patients deeper: paired Wilcoxon two-sided P = 2/1024 = 0.002 (lung, Fig. S22)."""
    near, pat, arm = _cohort(10, [16 + i for i in range(10)], [4] * 10)
    r = cv.within_patient_occupancy(near, pat, arm, "tumour", "normal", region="anchor")
    assert (r["n_patients"], r["n_deeper"], r["n_lower"]) == (10, 10, 0)
    assert r["wilcoxon_two_sided"] == pytest.approx(2 / 1024)
    assert r["sign_test_one_sided"] == pytest.approx(1 / 1024)
    assert r["median_normal"] == pytest.approx(0.1) and r["median_malignant"] > 0.4


def test_sign_test_conventions_match_the_paper():
    """The paper's sign tests are one-sided: 9 of 10 -> 0.011; 15 of 22 -> 0.07."""
    near, pat, arm = _cohort(10, [24 - 2 * i if i < 9 else 1 for i in range(10)], [4] * 9 + [10])
    r = cv.within_patient_occupancy(near, pat, arm, "tumour", "normal", region="anchor")
    assert (r["n_deeper"], r["n_lower"]) == (9, 1)
    assert r["sign_test_one_sided"] == pytest.approx(0.011, abs=0.001)
    mal = [10 + (i % 7) for i in range(15)] + [1] * 7
    nor = [4] * 15 + [12 + i % 5 for i in range(7)]
    near, pat, arm = _cohort(22, mal, nor)
    r = cv.within_patient_occupancy(near, pat, arm, "tumour", "normal", region="anchor")
    assert (r["n_deeper"], r["n_lower"]) == (15, 7)
    assert r["sign_test_one_sided"] == pytest.approx(0.07, abs=0.005)
    assert r["sign_test_two_sided"] == pytest.approx(2 * r["sign_test_one_sided"], rel=1e-6)


def test_within_patient_applies_the_20_cell_floor():
    near, pat, arm = _cohort(5, [16] * 5, [4] * 5, n_cells=40)
    small = pd.concat([near, pd.Series(["cytotrophoblasts"] * 10)], ignore_index=True)
    pat2 = pd.concat([pat, pd.Series(["P99"] * 10)], ignore_index=True)
    arm2 = pd.concat([arm, pd.Series(["tumour"] * 10)], ignore_index=True)
    r = cv.within_patient_occupancy(small, pat2, arm2, "tumour", "normal", region="anchor")
    assert r["n_patients"] == 5 and "P99" not in r["per_patient"].index
    r10 = cv.within_patient_occupancy(small, pat2, arm2, "tumour", "normal", region="anchor", min_cells=5)
    assert r10["n_patients"] == 5            # P99 has no normal cells at all


def test_within_patient_with_no_eligible_patient_is_empty_not_an_error():
    near, pat, arm = _cohort(3, [5] * 3, [5] * 3, n_cells=10)
    r = cv.within_patient_occupancy(near, pat, arm, "tumour", "normal", region="anchor")
    assert r["n_patients"] == 0 and np.isnan(r["wilcoxon_two_sided"])


def test_complexity_matching_removes_a_depth_driven_difference():
    """Control of Note S-F: resample the malignant arm to the normal arm's detected-gene distribution."""
    rng = np.random.default_rng(0)
    labels, patient, arm, genes = [], [], [], []
    for p in range(6):
        g_m = rng.integers(1000, 5000, 60)
        g_n = rng.integers(1000, 3000, 60)
        labels += ["cytotrophoblasts" if g > 3000 else "hepatocytes" for g in g_m]    # core only in deep cells
        labels += ["hepatocytes"] * 60
        patient += [f"P{p}"] * 120
        arm += ["tumour"] * 60 + ["normal"] * 60
        genes += list(g_m) + list(g_n)
    idx = pd.RangeIndex(len(labels))
    near, pat, arm_s, ng = pd.Series(labels, idx), pd.Series(patient, idx), pd.Series(arm, idx), pd.Series(genes, idx)
    raw = cv.within_patient_occupancy(near, pat, arm_s, "tumour", "normal", region="anchor")
    matched = cv.within_patient_occupancy(near, pat, arm_s, "tumour", "normal", region="anchor",
                                          n_genes=ng, match_complexity=True)
    assert raw["mean_delta"] > 0.3
    assert matched["matched_delta"] < 0.05
    with pytest.raises(ValueError):
        cv.within_patient_occupancy(near, pat, arm_s, "tumour", "normal", match_complexity=True)


# ── recurrence ranking ───────────────────────────────────────────────────────

def test_recurrence_ranking_counts_cohorts_where_a_cell_is_in_the_top_five():
    def cohort(extra):
        return pd.Series(["cytotrophoblasts"] * 50 + ["megakaryocytes"] * 30 + [extra] * 20)
    cohorts = {f"c{i}": cohort("hepatocytes") for i in range(4)}
    cohorts["c4"] = cohort("astrocytes")
    out = cv.recurrence_ranking(cohorts)
    assert out.loc["cytotrophoblasts", "recurrence"] == 5
    assert out.loc["hepatocytes", "recurrence"] == 4 and out.loc["astrocytes", "recurrence"] == 1
    assert out.loc["cytotrophoblasts", "occupancy_pct"] == pytest.approx(50.0)
    assert out.loc["cytotrophoblasts", "median_depth_pct"] == pytest.approx(50.0)
    assert out["occupancy_pct"].sum() == pytest.approx(100.0)
    assert "erythrocytes" not in out.index and len(out) == 152
    assert out.index[0] == "cytotrophoblasts"


# ── Table S7 check and module contributions ──────────────────────────────────

def test_pseudobulk_calibration_does_not_reproduce_table_s7():
    """The paper's Table S7 is computed on the single-cell-scaled calibration; the bundled
    pseudobulk reference is correlated with it but does not match it."""
    res = cv.check_calibration_against_table_s7("pseudobulk")
    assert set(res) == {"megakaryocyte_r", "cytotrophoblast_r", "max_abs_diff", "mean_abs_diff",
                        "max_abs_diff_checked", "n_checked", "n_exact", "known_deviations", "matches"}
    assert 0.5 < res["megakaryocyte_r"] < 0.95 and 0.5 < res["cytotrophoblast_r"] < 0.95
    assert res["matches"] is False and res["n_exact"] < 10 and res["n_checked"] == 86


def test_a_calibration_built_from_table_s7_values_would_match():
    cal = pct.load_calibration("pseudobulk")
    # replace the two core reference rows so that their z-scores equal Table S7 exactly
    ref = cal.reference_scores.copy()
    z = cal.standardize(ref)
    z.loc["megakaryocytes", S7.index] = S7["megakaryocyte_z"].values
    z.loc["cytotrophoblasts", S7.index] = S7["cytotrophoblast_z"].values
    ref2 = z * cal.std + cal.mean
    cal2 = pct.Calibration(cal.name, cal.module_order, cal.mean.values, cal.std.values, cal.loadings, reference_scores=ref2)
    res = cv.check_calibration_against_table_s7(cal2)
    assert res["matches"] is True and res["n_exact"] == 86 and res["known_deviations"]["cGAS-STING"] == pytest.approx((0.0, 0.0), abs=1e-9)
    # the known deviation is left out of the verdict but never hidden
    assert cv.check_calibration_against_table_s7(cal, known_deviations=())["n_checked"] == 88


def test_module_contributions_sum_to_the_cosine_with_the_signature():
    cal = pct.load_calibration("pseudobulk")
    sig = cv.core_signature(cal)
    assert np.linalg.norm(sig.values) == pytest.approx(1.0)
    row = cal.reference_scores.loc["gastric progenitor cells"]
    contrib = cv.module_contributions(row, cal)
    z = cal.standardize(row.to_frame().T).iloc[0]
    expected = float(z @ sig / np.linalg.norm(z))
    assert contrib.sum() == pytest.approx(expected)
    # the core itself is close to its own signature
    assert cv.module_contributions(cal.reference_scores.loc["cytotrophoblasts"], cal).sum() > 0.7
    prof = cv.core_profile(cal)
    assert list(prof.columns) == ["cytotrophoblasts", "megakaryocytes", "mean"] and len(prof) == 44
