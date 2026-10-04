"""Synthetic checks of scripts/19_gate0prime_coretention.py.

Worlds (genes gained at an origin age and then lost along lineages):
  tautology   modules = genes of similar origin, independent losses. The old P5-style criterion
              (intact-species count, ancient > young) passes; gate 0' must NOT pass.
  two_layer   ancient modules lose their genes together; young modules are assembled from genes whose
              losses fall in complementary species. Gate 0' must give PASS_TWO_LAYERS.
  ancient     only the ancient co-loss -> PASS_ANCIENT_UNITS_ONLY.
"""

import importlib.util
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DIV = np.array([90, 90, 180, 310, 435, 435, 684, 684, 797, 797, 950, 1105, 1500, 1500, 1500], float)
ORIGINS = [1500, 1105, 797, 684, 435]


def _load():
    spec = importlib.util.spec_from_file_location("g0", ROOT / "scripts" / "19_gate0prime_coretention.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.N_NULL = 1000
    return m


def _world(seed, kind, n_mod=12, size=6, n_bg=2000, loss=0.12):
    rng = np.random.default_rng(seed)
    rows, mods_cat = [], {}

    def gene(origin, lost=()):
        p = (DIV <= origin).astype(np.int8)
        p[(rng.random(len(DIV)) < loss) & (p == 1)] = 0
        p[list(lost)] = 0
        return p

    for i in range(n_bg):
        rows.append(gene(ORIGINS[i % len(ORIGINS)]))
    k = 0
    for origin in ORIGINS:
        for j in range(n_mod):
            name = f"M{origin}_{j}"
            genes = []
            present = np.where(DIV <= origin)[0]
            if kind in ("two_layer", "ancient") and origin >= 1500:
                shared = present[rng.random(len(present)) < 0.35]          # losses of the whole module
                for _ in range(size):
                    p = (DIV <= origin).astype(np.int8)
                    p[shared] = 0
                    rows.append(p)
                    genes.append(len(rows) - 1)
            elif kind == "two_layer" and origin <= 797:
                # complementary: each gene lost in its own slice of the species where the module exists
                order = rng.permutation(present)
                slices = np.array_split(order, size)
                for sl in slices:
                    p = (DIV <= origin).astype(np.int8)
                    p[sl] = 0
                    rows.append(p)
                    genes.append(len(rows) - 1)
            else:
                for _ in range(size):
                    rows.append(gene(origin))
                    genes.append(len(rows) - 1)
            mods_cat[name] = genes
            k += 1
    P = np.vstack(rows)
    names = np.array([f"g{i}" for i in range(len(P))])
    mods = {m: [names[i] for i in g] for m, g in mods_cat.items()}
    return names, P, mods


def _gate(g0, seed, kind):
    per = {}
    for c in ("C0", "C2", "C3", "C4"):
        names, P, mods = _world(seed + hash(c) % 97, kind)
        rng = np.random.default_rng(1)
        rows = g0.analyse_catalog(mods, names, P, DIV, P.sum(1), rng)
        per[c] = g0.criteria(rows, g0.calibration(mods, names, P, DIV, P.sum(1), rng))
        per[c]["_rows"] = rows
    return g0.verdict(per, ["C2", "C3", "C4"]), per


def _old_g01(rows):
    age = np.array([r["age_mya"] for r in rows])
    n = np.array([r["observed_intact"] for r in rows])
    return stats.mannwhitneyu(n[age >= 1500], n[age <= 800], alternative="greater").pvalue


def test_tautology_world_old_criterion_passes_new_does_not():
    g0 = _load()
    v, per = _gate(g0, 11, "tautology")
    assert all(_old_g01(per[c]["_rows"]) < 0.05 for c in ("C2", "C3", "C4"))
    assert v in ("FAIL", "NOT_EVALUABLE"), (v, {c: per[c].get("p_L1_ancient_coherent") for c in per})


def test_two_layer_world_passes():
    g0 = _load()
    v, per = _gate(g0, 21, "two_layer")
    assert v == "PASS_TWO_LAYERS", (v, {c: (per[c]["L1"], per[c]["L2"], per[c]["L3"]) for c in per})


def test_ancient_only_world():
    g0 = _load()
    v, per = _gate(g0, 31, "ancient")
    assert v == "PASS_ANCIENT_UNITS_ONLY", (v, {c: (per[c]["L1"], per[c]["L2"], per[c]["L3"]) for c in per})


def test_tautology_false_positive_rate():
    g0 = _load()
    l2 = 0
    for s in range(5):
        v, per = _gate(g0, 100 + 13 * s, "tautology")
        assert not v.startswith("PASS")
        l2 += sum(bool(per[c]["L2"]) for c in ("C2", "C3", "C4"))
    assert l2 <= 2                      # 15 catalog-level L2 tests under no effect


def test_calibration_flags_a_biased_null():
    g0 = _load()
    names, P, mods = _world(7, "tautology", n_bg=300)    # small pool -> widened classes
    rng = np.random.default_rng(2)
    rows = g0.analyse_catalog(mods, names, P, DIV, P.sum(1), rng)
    c = g0.criteria(rows, g0.calibration(mods, names, P, DIV, P.sum(1), rng))
    assert c["widened_frac"] > 0 and "calib_mean_z" in c
