"""Synthetic checks of the locked loop-integrity analysis (scripts/14_...).

Planted truth:
  - positive world: break KDs raise output dispersion (and the mean), input_down KDs lower the mean;
    the analysis must reach POSITIVE.
  - mean-only world: break KDs raise the mean only, no extra dispersion; the analysis must NOT reach POSITIVE.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ad = pytest.importorskip("anndata")
import perceptome as pct  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("li", ROOT / "scripts" / "14_loop_integrity_prereg_v1.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.N_NULL, m.N_BOOT, m.N_PERM, m.N_PSEUDO, m.PSEUDO_CELLS = 200, 60, 2000, 12, 100
    return m


def _world(tmp_path, li, dispersion_effect):
    rng = np.random.default_rng(1)
    loops = pct.load_loops()["loops"]
    loop_genes = sorted({g for e in loops.values() for k in ("feedback_core", "loop_targets") for g in e[k]})
    genes = loop_genes + [f"POOL{i}" for i in range(500)]
    gi = {g: i for i, g in enumerate(genes)}
    base_mu = rng.lognormal(mean=0.0, sigma=1.0, size=len(genes))
    base_phi = np.full(len(genes), 0.2)

    aset = li.analysis_set(li.read_contrasts())
    labels, n_kd = [], 90
    kd_effect = {}
    for m in li.PRIMARY:
        out = [gi[g] for g in pct.loop_output_genes(m) if g in gi]
        for g in aset.get(m, {}).get("break", []):
            kd_effect[g] = (out, 1.4, 6.0 if dispersion_effect else 1.0)
        for g in aset.get(m, {}).get("input_down", []):
            kd_effect[g] = (out, 0.5, 1.0)
    kds = sorted({u for m in aset for c in aset[m].values() for u in c} | {"MBTPS1", "MBTPS2", "ATF6", "SCAP"})
    for g in kds:
        labels += [g] * n_kd
    labels += ["non-targeting"] * 2500
    n = len(labels)
    gem = rng.choice(["g1", "g2", "g3", "g4"], size=n)
    size = rng.lognormal(0, 0.25, size=n) * 6.0
    X = np.empty((n, len(genes)), dtype=np.float32)
    for c in range(n):
        mu, phi = base_mu.copy(), base_phi.copy()
        eff = kd_effect.get(labels[c])
        if eff:
            mu[eff[0]] *= eff[1]
            phi[eff[0]] *= eff[2]
        lam = rng.gamma(1 / phi, mu * size[c] * phi)
        X[c] = rng.poisson(lam)
    import pandas as pd
    obs = pd.DataFrame({"gene": labels, "gem_group": gem}, index=[str(i) for i in range(n)])
    var = pd.DataFrame({"gene_name": genes}, index=[f"ENSG{i}" for i in range(len(genes))])
    path = tmp_path / "fake.h5ad"
    ad.AnnData(X, obs=obs, var=var).write_h5ad(path)
    return path


def test_positive_world_detected(tmp_path):
    li = _load()
    v = li.run(_world(tmp_path, li, dispersion_effect=True), tmp_path / "out")
    assert v["G2"]["pass"], v["G2"]
    assert len(v["loops_passing_G1"]) >= 3, v["G1"]
    assert v["decision"].startswith("POSITIVE"), v


def test_mean_only_world_not_positive(tmp_path):
    li = _load()
    v = li.run(_world(tmp_path, li, dispersion_effect=False), tmp_path / "out")
    assert not v["decision"].startswith("POSITIVE"), v
