"""Synthetic checks of the locked homeostatic-locality analysis (scripts/18_...).

Worlds (two cell lines each):
  - effect:      knockdowns in interoceptive modules move their own targets more -> CONFIRMED
  - null:        every module equally local                                      -> not CONFIRMED / PARTIAL
  - detectability confound: interoceptive targets are highly expressed and z-noise grows
    with expression, with no true difference -> the matched null must keep it from CONFIRMED
"""

import importlib.util
from pathlib import Path

import numpy as np

import perceptome as pct

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("hl", ROOT / "scripts" / "18_homeostatic_locality_prereg_v2.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.N_NULL = 100
    return m


def _world(hl, seed, delta, confound=False, base_effect=1.0):
    rng = np.random.default_rng(seed)
    mods = pct.list_modules()
    act = {m: set(pct.get_genes(m, "activity")) for m in mods}
    intero = {m for m in mods if hl.category(m) == hl.GROUP}
    cat_genes = sorted(set().union(*act.values()) | set().union(*(hl.machinery(m) for m in mods)))
    genes = np.array(cat_genes + [f"POOL{i}" for i in range(3000)])
    gi = {g: i for i, g in enumerate(genes)}
    umi = rng.lognormal(-0.5, 1.2, len(genes))
    if confound:
        for m in intero:
            for g in act[m]:
                umi[gi[g]] = rng.lognormal(2.5, 0.3)
    noise_sd = 1.0 + (0.6 * np.log10(umi + 1) * 3 if confound else 0.0)
    kds = sorted(set().union(*(hl.machinery(m) for m in mods)))
    Z = np.zeros((len(kds), len(genes)), dtype=np.float32)
    for k, kd in enumerate(kds):
        z = rng.normal(0, 0.15, len(genes)) * noise_sd
        for m in mods:
            if kd in hl.machinery(m):
                e = base_effect + (delta if m in intero else 0.0)
                for g in act[m]:
                    z[gi[g]] += rng.normal(0, 0.15 * e) * (noise_sd[gi[g]] if confound else 1.0)
        Z[k] = z
    return {"Z": Z, "genes": genes, "kd": np.array(kds), "n_cells": np.full(len(kds), 120),
            "fold": np.full(len(kds), 0.3), "kd_ctrl_mean_umi": np.full(len(kds), 1.0),
            "ctrl_mean_umi": umi, "ctrl_mean_log": np.log1p(umi)}


def _run(hl, seed, **kw):
    res = {ln: hl.analyse_line(_world(hl, seed + i, **kw))[0]
           for i, ln in enumerate((hl.PRIMARY_LINE, hl.REPLICATION_LINE))}
    return hl.verdict(res), res


def test_effect_world_confirms():
    hl = _load()
    v, res = _run(hl, 10, delta=3.0)
    assert res[hl.PRIMARY_LINE]["primary"]["evaluable"]
    assert v["verdict"] == "CONFIRMED", v


def test_null_world_does_not_confirm():
    hl = _load()
    hits = 0
    for s in range(6):
        v, _ = _run(hl, 100 + 10 * s, delta=0.0)
        hits += v["verdict"].startswith("CONFIRMED")
    assert hits == 0


def test_detectability_confound_does_not_confirm():
    hl = _load()
    v, res = _run(hl, 200, delta=0.0, confound=True)
    assert not v["verdict"].startswith("CONFIRMED"), (v, res[hl.PRIMARY_LINE]["primary"])


def test_off_module_control_is_quiet_in_effect_world():
    hl = _load()
    _, res = _run(hl, 300, delta=3.0)
    assert res[hl.PRIMARY_LINE]["N1_off_module"].get("p_one_sided", 1) >= hl.ALPHA


def test_builder_on_small_h5ad(tmp_path):
    """End-to-end: CSR h5ad -> pseudobulk; z and fold match a direct computation."""
    import argparse

    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    spec = importlib.util.spec_from_file_location("pb", ROOT / "scripts" / "17_xatlas_pseudobulk.py")
    pb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pb)
    rng = np.random.default_rng(0)
    genes = ["HMGCR", "SCAP", "FASN", "LDLR"] + [f"G{i}" for i in range(40)]
    n_ctrl, n_kd = 600, 120
    lab = ["NTC"] * n_ctrl + ["SCAP"] * n_kd + ["HMGCR"] * n_kd + ["OTHER"] * 50
    mu = rng.uniform(0.5, 5, len(genes))
    C = rng.poisson(mu, (len(lab), len(genes))).astype(np.float32)
    C[n_ctrl:n_ctrl + n_kd, 1] = rng.poisson(mu[1] * 0.2, n_kd)        # SCAP knocked down to ~20 %
    C[n_ctrl:n_ctrl + n_kd, 3] = rng.poisson(mu[3] * 3, n_kd)          # LDLR responds
    perm = rng.permutation(len(lab))                                   # cells not sorted by target
    C, lab = C[perm], np.array(lab)[perm]
    a = ad.AnnData(X=sp.csr_matrix(C), obs=pd.DataFrame({"gene_target": pd.Categorical(lab),
                                                          "batch": pd.Categorical(["b1"] * len(lab))},
                                                         index=[f"c{i}" for i in range(len(lab))]),
                   var=pd.DataFrame(index=genes))
    f = tmp_path / "toy.h5ad"
    a.write_h5ad(f)
    out = tmp_path / "pb"
    pb.build(argparse.Namespace(h5ad=str(f), out=str(out), target_col="gene_target", control_label="NTC",
                                batch_col="batch", gene_col=None, filter_col=None))
    d = np.load(out / "pseudobulk.npz")
    kd = list(d["kd"].astype(str))
    assert set(kd) == {"SCAP", "HMGCR"}
    i = kd.index("SCAP")
    assert d["n_cells"][i] == n_kd
    assert 0.1 < d["fold"][i] < 0.35
    tot = C.sum(1, keepdims=True)
    L = np.log1p(C / tot * 1e4)
    ctrl, k = lab == "NTC", lab == "SCAP"
    z = ((L[k] - L[ctrl].mean(0)) / np.maximum(L[ctrl].std(0), pb.SIGMA_FLOOR)).mean(0)
    assert np.allclose(d["Z"][i], z, atol=1e-3)
    assert d["Z"][i][3] > 1.0


def test_extractor_subset_gives_identical_pseudobulk(tmp_path):
    """Amendment 1: script 17 on the 17a subset == script 17 on the full file (incl. control sampling)."""
    import argparse

    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    def load(name, fname):
        spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / fname)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    pb, ex = load("pbA", "17_xatlas_pseudobulk.py"), load("exA", "17a_extract_subset.py")
    pb.N_CTRL = ex.pb.N_CTRL = 150                       # force real sampling from 400 controls
    ex.ROWS_PER_CHUNK = 37
    rng = np.random.default_rng(3)
    genes = ["HMGCR", "SCAP", "FASN", "LDLR"] + [f"G{i}" for i in range(30)]
    lab = np.array(["NTC"] * 400 + ["SCAP"] * 90 + ["HMGCR"] * 70 + ["OTHER"] * 300)[rng.permutation(860)]
    C = rng.poisson(rng.uniform(0.3, 4, len(genes)), (len(lab), len(genes))).astype(np.float32)
    filt = rng.random(len(lab)) > 0.1
    obs = pd.DataFrame({"gene_target": pd.Categorical(lab),
                        "sample": pd.Categorical(rng.choice(["s1", "s2"], len(lab))),
                        "pass_guide_filter": filt}, index=[f"c{i}" for i in range(len(lab))])
    full = tmp_path / "full.h5ad"
    ad.AnnData(X=sp.csr_matrix(C), obs=obs, var=pd.DataFrame(index=genes)).write_h5ad(full)
    flags = dict(target_col="gene_target", control_label="NTC", batch_col="sample", gene_col=None,
                 filter_col="pass_guide_filter")
    sub = tmp_path / "sub.h5ad"
    import sys
    argv = sys.argv
    sys.argv = ["17a", "--h5ad", str(full), "--out", str(sub), "--work", str(tmp_path / "w"),
                "--target_col", "gene_target", "--control_label", "NTC", "--batch_col", "sample",
                "--filter_col", "pass_guide_filter"]
    try:
        ex.main()
    finally:
        sys.argv = argv
    pb.MIN_CTRL_PER_BATCH = 10
    pb.build(argparse.Namespace(h5ad=str(full), out=str(tmp_path / "pf"), **flags))
    pb.build(argparse.Namespace(h5ad=str(sub), out=str(tmp_path / "ps"), **flags))
    a, b = np.load(tmp_path / "pf" / "pseudobulk.npz"), np.load(tmp_path / "ps" / "pseudobulk.npz")
    for k in ("kd", "n_cells", "genes"):
        assert np.array_equal(a[k], b[k]), k
    for k in ("Z", "fold", "ctrl_mean_umi", "ctrl_mean_log"):
        assert np.allclose(a[k], b[k], equal_nan=True), k
