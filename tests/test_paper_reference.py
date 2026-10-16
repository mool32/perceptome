"""Statements of the paper (Spiro 2026) that the bundled reference reproduces.

Everything here is checked against the *bundled pseudobulk reference* (reference_v03).
Statements computed on the single-cell-scaled calibration (Fig. 1B variance labels,
Fig. S8, Table S7 z-scores) are not reproducible from it; they are listed in
docs/PAPER_MAP.md instead of being asserted here.
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
