"""Statements of the paper (Spiro 2026) that the bundled reference reproduces.

The first part is checked against the *bundled pseudobulk reference* (reference_v03): the eigenspace
construction (Kaiser 9, 73.8%, stability, gene lists). The last part is checked against the *bundled
single-cell-scaled calibration*, the reference the paper re-derives for per-cell data (Fig. 1B shares,
Fig. S8, Table S7, the convergent neighbourhood, the medoid, Note S-N iii). docs/PAPER_MAP.md lists what
is and is not reproduced.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import perceptome as pct
from perceptome.eigenspace.project import _load_reference

_DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def ref():
    return _load_reference()


@pytest.fixture(scope="module")
def table_s9():
    return json.loads((_DATA / "paper_table_s9_core_genes.json").read_text(encoding="utf-8"))["modules"]


# ── Methods / Note S1: eigenspace construction ───────────────────────────────

def test_nine_components_retained_by_the_kaiser_criterion(ref):
    assert ref["n_kaiser"] == ref["n_pcs"] == 9
    ev = np.array(ref["eigenvalues"])
    assert (ev[:9] > 1).all() and ev[9] <= 1


def test_nine_components_explain_73_8_percent_of_variance(ref):
    """Note S1: 'Nine components were retained by the Kaiser criterion (73.8% of variance)'."""
    cum = sum(ref["var_explained"][:9])
    assert 0.738 <= cum < 0.739


def test_bootstrap_stability_tiers_match_fig_s1(ref):
    """Fig. S1: PC1-PC6 stable (0.63-0.99), PC7-PC9 probable (0.40-0.54); PC8 dips below PC9."""
    s = [ref["stability"][f"PC{i}"] for i in range(1, 10)]
    assert (round(min(s[:6]), 2), round(max(s[:6]), 2)) == (0.63, 0.99)
    assert (round(min(s[6:]), 2), round(max(s[6:]), 2)) == (0.40, 0.54)
    assert s[7] < s[8]
    assert [ref["pc_confidence"][f"PC{i}"] for i in range(1, 10)] == ["stable"] * 6 + ["probable"] * 3


def test_reference_is_reproduced_by_rebuild_from_the_bundled_matrix(ref):
    R = pct.load_hpa_perceptivity()["R"]
    out = pct.rebuild(R, bootstrap_n=0)
    np.testing.assert_allclose(out["eigenvalues"], ref["eigenvalues"], atol=1e-8)
    for j in range(9):
        a = np.array([out["loadings"][f"PC{j + 1}"][m] for m in ref["module_order"]])
        b = np.array([ref["loadings"][f"PC{j + 1}"][m] for m in ref["module_order"]])
        assert abs(a @ b) == pytest.approx(1.0, abs=1e-8)


def test_npas4_leave_one_out_reproduces_all_nine_components(ref):
    """Note S2: rebuilding on 43 modules gives per-PC eigenvector cosine >= 0.95 (mean 0.99 over PC1-PC6)."""
    R = pct.load_hpa_perceptivity()["R"]
    mods = [m for m in ref["module_order"] if m != "NPAS4"]
    out = pct.rebuild(R[mods], bootstrap_n=0)
    cos = []
    for j in range(9):
        a = np.array([out["loadings"][f"PC{j + 1}"][m] for m in mods])
        b = np.array([ref["loadings"][f"PC{j + 1}"][m] for m in mods])
        cos.append(abs(a @ b) / np.linalg.norm(b))
    assert min(cos) >= 0.95
    assert np.mean(cos[:6]) == pytest.approx(0.99, abs=0.01)


def test_reference_has_154_cell_types_and_44_modules():
    ref = pct.load_hpa_perceptivity()
    assert ref["R"].shape == (154, 44)


# ── Table S9: the catalog is the explicit gene list of the paper ─────────────

def test_catalog_core_genes_are_table_s9(table_s9):
    cat = pct.load_catalog()["modules"]
    assert set(cat) == set(table_s9)
    for m, genes in table_s9.items():
        assert cat[m]["core_genes"] == genes, m


def test_catalog_gene_counts_match_the_paper(table_s9):
    """219 gene slots reduce to 206 unique genes; Calcium is the only 4-gene module (Table S9)."""
    slots = [g for genes in table_s9.values() for g in genes]
    assert len(slots) == 219 and len(set(slots)) == 206
    sizes = {m: len(g) for m, g in table_s9.items()}
    assert sizes["Calcium"] == 4
    assert sum(1 for n in sizes.values() if n == 5) == 43


def test_twelve_genes_are_shared_across_17_modules(table_s9):
    """Note S-L (ii): twelve genes appear in more than one module, touching 17 of the 44 modules."""
    holders = {}
    for m, genes in table_s9.items():
        for g in genes:
            holders.setdefault(g, []).append(m)
    shared = {g: ms for g, ms in holders.items() if len(ms) > 1}
    assert set(shared) == {"ADIPOQ", "CALM1", "CAMK2A", "CCND1", "DUSP1", "FKBP5",
                           "FOXO1", "GSK3B", "PGR", "SGK1", "STAT1", "XBP1"}
    assert len({m for ms in shared.values() for m in ms}) == 17


def test_definitional_module_pairs_share_a_gene(table_s9):
    """Note S-L (ii): UPR-ATF6/UPR-IRE1 share XBP1 and JAK-STAT/Type I IFN share STAT1."""
    assert "XBP1" in table_s9["UPR-ATF6"] and "XBP1" in table_s9["UPR-IRE1"]
    assert "STAT1" in table_s9["JAK-STAT"] and "STAT1" in table_s9["Type I IFN"]
    assert not set(table_s9["ERK/MAPK"]) & set(table_s9["NF-κB"])


# ── reproduction report (script 15) ──────────────────────────────────────────

def _load_report_script():
    import importlib.util
    path = Path(__file__).resolve().parents[1] / "scripts" / "15_paper_reproduction_report.py"
    spec = importlib.util.spec_from_file_location("paper_report", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_reproduction_report_on_the_pseudobulk_calibration():
    """The pseudobulk reference reproduces the eigenspace construction (Kaiser 9, 73.8%) and nothing that
    belongs to the single-cell-scaled reference (Fig. S8, Table S7, Fig. 1B shares, neighbourhood, medoid)."""
    rows = _load_report_script().rows_for(pct.load_calibration("pseudobulk"))
    verdict = {r[1]: r[4] for r in rows}
    assert len(rows) == 12
    assert [v for v in verdict.values()].count("match") == 3
    assert verdict["components with eigenvalue > 1 (Kaiser)"] == "match"
    assert verdict["variance of the nine components"] == "match"
    assert verdict["random pairs: mean cosine"] == "match"
    assert verdict["module pairs with |r| > 0.7 of 946 (Fig. S8)"] == "differs"
    assert verdict["four-cell neighbourhood: mutual cosine"] == "differs"


def test_reproduction_report_on_the_single_cell_scaled_calibration():
    """10 of 12 statements are reproduced; the two that are not are the cGAS-STING row of Table S7 (flagged by
    the paper's Note S-Q) and the 0.99 mean distance among the origins (1.000 here)."""
    rows = _load_report_script().rows_for(pct.load_calibration("single_cell_scaled"))
    verdict = {r[1]: r[4] for r in rows}
    assert len(rows) == 12
    assert sorted(k for k, v in verdict.items() if v == "differs") == [
        "Table S7, cGAS-STING row (flagged in Note S-Q)", "mean pairwise cosine distance among the 25 origins"]
    assert verdict["medoid of the normal manifold (cosine distance, 9-D)"] == "match"


# ── the single-cell-scaled calibration (the paper's reference for per-cell data) ──────────────────────

PAPER_PAIRS = {frozenset(p) for p in [("UPR-ATF6", "UPR-IRE1"), ("JAK-STAT", "Type I IFN"), ("ERK/MAPK", "NF-κB"),
                                      ("UPR-PERK", "p53"), ("NF-κB", "UPR-PERK"), ("BMP", "Hippo")]}


@pytest.fixture(scope="module")
def sc():
    return pct.load_calibration("single_cell_scaled")


@pytest.fixture(scope="module")
def s7():
    from perceptome import convergence as cv
    return cv.load_convergence_reference()["tables"]["S7"]


def test_the_single_cell_calibration_is_the_papers_construction(sc):
    """Methods: the reference is re-derived under single_cell_scaled; the 44 x 9 loadings stay fixed."""
    ps = pct.load_calibration("pseudobulk")
    assert sc.name == "single_cell_scaled" and sc.ddof == 1 and sc.n_pcs == 9
    assert list(sc.reference_scores.index) == list(ps.reference_scores.index) and len(sc.reference_scores) == 154
    np.testing.assert_allclose(sc.loadings.values, ps.loadings.values)
    np.testing.assert_allclose(sc.std.values, sc.reference_scores[sc.module_order].std(ddof=1).values)
    # per-cell scores sit on a different scale from the population averages: the reason it is needed
    assert (ps.reference_scores.mean().values > sc.reference_scores.mean().values).all()


def test_table_s7_core_z_scores_are_reproduced_to_the_printed_digit_for_43_modules(sc, s7):
    from perceptome import convergence as cv
    res = cv.check_calibration_against_table_s7(sc)
    assert res["matches"] and res["n_checked"] == 86
    assert res["n_exact"] >= 85 and res["max_abs_diff_checked"] < 0.006
    assert res["megakaryocyte_r"] > 0.999 and res["cytotrophoblast_r"] > 0.999
    # reported, not hidden: the paper's own Note S-Q flags the cGAS-STING module (three genes, no sensor)
    assert set(res["known_deviations"]) == set(cv.S7_KNOWN_DEVIATIONS) == {"cGAS-STING"}
    dm, dc = res["known_deviations"]["cGAS-STING"]
    assert dm == pytest.approx(-0.110, abs=0.005) and dc == pytest.approx(-0.135, abs=0.005)
    assert res["max_abs_diff"] == pytest.approx(0.135, abs=0.005)


def test_table_s7_target_distinctiveness_is_the_z_of_the_core_against_the_other_152_cells(sc, s7):
    from perceptome import convergence as cv
    d = cv.target_distinctiveness(sc).reindex(s7.index)
    diff = (d - s7["target_distinctiveness"]).abs().drop("cGAS-STING")
    assert diff.max() < 0.003 and int((diff < 0.0006).sum()) >= 41
    np.testing.assert_allclose(s7["core_suppression"], -s7["target_distinctiveness"])
    assert d.idxmax() == "cGAS-STING" and d.idxmin() == "Circadian"                # the one inverted sign
    assert float(s7["cancer_contribution"].sum()) == pytest.approx(0.711, abs=0.0006)


def test_fig_s8_module_correlations_come_from_the_single_cell_matrix(sc):
    """Fig. S8: six of 946 module pairs exceed |r| > 0.7; mean |r| 0.21. The pseudobulk matrix gives 20 and 0.27."""
    def stats(cal):
        X = cal.reference_scores[cal.module_order].values
        r = np.corrcoef(X.T)
        iu = np.triu_indices(len(r), 1)
        big = {frozenset((cal.module_order[i], cal.module_order[j])) for i, j in zip(*iu) if abs(r[i, j]) > 0.7}
        return big, float(np.abs(r[iu]).mean())
    big, mean_abs = stats(sc)
    assert len(big) == 6 and big == PAPER_PAIRS
    assert mean_abs == pytest.approx(0.21, abs=0.005)
    big_ps, mean_ps = stats(pct.load_calibration("pseudobulk"))
    assert len(big_ps) == 20 and mean_ps == pytest.approx(0.27, abs=0.01)


def test_fig_1b_variance_labels_are_axis_shares_of_the_single_cell_projection(sc):
    """'the first component (31% of variance) ... the second (21%)': the share of each fixed axis among the
    nine retained ones when the 154 cell types are projected under the single-cell calibration."""
    v = sc.centroids.var(axis=0, ddof=0)
    share = 100 * v / v.sum()
    assert share.iloc[0] == pytest.approx(31, abs=0.5) and share.iloc[1] == pytest.approx(21, abs=0.5)
    vp = pct.load_calibration("pseudobulk").centroids.var(axis=0, ddof=0)
    assert 100 * vp.iloc[0] / vp.sum() > 38                       # not the pseudobulk reading (39.5%)


def test_the_convergent_neighbourhood_has_mutual_cosine_0_77_against_minus_0_01_for_random_pairs(sc):
    import itertools
    from perceptome import convergence as cv
    names = list(sc.centroids.index)
    U = sc.centroids.values / np.linalg.norm(sc.centroids.values, axis=1, keepdims=True)
    S = U @ U.T
    nb = cv.region_members("neighbourhood")
    mutual = np.mean([S[names.index(a), names.index(b)] for a, b in itertools.combinations(nb, 2)])
    iu = np.triu_indices(len(S), 1)
    assert mutual == pytest.approx(0.77, abs=0.01)
    assert S[iu].mean() == pytest.approx(-0.01, abs=0.01)
    Up = pct.load_calibration("pseudobulk").centroids.values
    Up = Up / np.linalg.norm(Up, axis=1, keepdims=True)
    namesp = list(pct.load_calibration("pseudobulk").centroids.index)
    mutual_p = np.mean([(Up @ Up.T)[namesp.index(a), namesp.index(b)] for a, b in itertools.combinations(nb, 2)])
    assert mutual_p == pytest.approx(0.52, abs=0.02)


def test_the_medoid_of_the_normal_manifold_is_the_alveolar_type_1_pneumocyte(sc):
    """'the true medoid — an alveolar type-1 pneumocyte': the cell type with the smallest summed cosine distance."""
    U = sc.centroids.values / np.linalg.norm(sc.centroids.values, axis=1, keepdims=True)
    names = list(sc.centroids.index)
    D = 1 - U @ U.T
    order = np.argsort(D.sum(axis=1))
    assert names[order[0]] == "alveolar cells type 1"
    assert names[order[1]] == "transitional alveolar cells"


def test_note_s_n_iii_the_rebuilt_space_differs_from_the_shipped_reference_by_r_0_964(sc):
    """Note S-N (iii): the eigenspace rebuilt from the normal cells and the shipped nine-dimensional reference are
    'near-identical but not identical (cell-to-cell geometry r = 0.964)'. On the single-cell matrix: the
    correlation of pairwise cosine distances between the rebuilt space and the shipped loadings."""
    from scipy.spatial.distance import pdist
    X = sc.reference_scores[sc.module_order].astype(float)
    Z = (X - X.mean()) / X.std(ddof=0)
    w, V = np.linalg.eigh(np.corrcoef(X.values.T))
    rebuilt = Z.values @ V[:, ::-1][:, :9]
    r = np.corrcoef(pdist(rebuilt, "cosine"), pdist(sc.centroids.values, "cosine"))[0, 1]
    assert r == pytest.approx(0.964, abs=0.003)
    ps = pct.load_calibration("pseudobulk")
    Xp = ps.reference_scores[ps.module_order].astype(float)
    Zp = (Xp - Xp.mean()) / Xp.std(ddof=0)
    wp, Vp = np.linalg.eigh(np.corrcoef(Xp.values.T))
    rp = np.corrcoef(pdist(Zp.values @ Vp[:, ::-1][:, :9], "cosine"), pdist(ps.centroids.values, "cosine"))[0, 1]
    assert rp == pytest.approx(1.0, abs=1e-6)                      # the pseudobulk rebuild is the shipped reference


def test_the_origins_are_near_orthogonal_under_the_single_cell_calibration(sc):
    """0.99 printed; 1.000 here (0.990 in the shape-only variant) — recorded as not reproduced in PAPER_MAP."""
    import itertools
    from perceptome import convergence as cv
    origins = list(cv.load_convergence_reference()["tables"]["S8"]["cell_of_origin"])
    names = list(sc.centroids.index)
    U = sc.centroids.values / np.linalg.norm(sc.centroids.values, axis=1, keepdims=True)
    S = U @ U.T
    d = np.mean([1 - S[names.index(a), names.index(b)] for a, b in itertools.combinations(origins, 2)])
    assert d == pytest.approx(1.000, abs=0.002)
    so = sc.shape_only_calibration()
    names_so = list(so.centroids.index)
    Uo = so.centroids.values / np.linalg.norm(so.centroids.values, axis=1, keepdims=True)
    So = Uo @ Uo.T
    d_so = np.mean([1 - So[names_so.index(a), names_so.index(b)] for a, b in itertools.combinations(origins, 2)])
    assert d_so == pytest.approx(0.990, abs=0.002)
