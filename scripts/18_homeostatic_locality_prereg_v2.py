"""Homeostatic-locality test, pre-registration v2 (docs/PREREG_homeostatic_locality_v2.md).

This script IS the pre-registered analysis, locked with the pre-registration
(sha256 in docs/PREREG_LOCK_v2.txt). Do not edit after the lock; any change is a
deviation and goes to docs/PREREG_homeostatic_locality_v2_DEVIATIONS.md.

Hypothesis: knocking down machinery of an interoceptive (internal-state-sensing) module
moves that module's own target genes more, relative to the rest of the transcriptome,
than knocking down machinery of a module that relays external signals.

Input: pseudobulk.npz written by scripts/17_xatlas_pseudobulk.py for each cell line.

Usage:
  python scripts/18_homeostatic_locality_prereg_v2.py \
      --line HCT116=results/xatlas_pb/HCT116 --line HEK293T=results/xatlas_pb/HEK293T \
      --out results/homeostatic_locality_v2
"""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perceptome as pct  # noqa: E402

# ---- locked parameters (§4-§7) ------------------------------------------------
MIN_CELLS = 50
MAX_FOLD = 0.50                  # target CP10k in knockdown cells / controls
SENS_MAX_FOLD = 0.60             # sensitivity only
MIN_KD_CTRL_UMI = 0.10           # target must be detectable in controls to judge the fold
MIN_GENE_CTRL_UMI = 0.01         # readout universe
MIN_ACTIVITY_GENES = 3
N_NULL = 500
N_DECILES = 10
ALPHA = 0.05
MIN_PER_GROUP = 8
PRIMARY_LINE = "HCT116"
REPLICATION_LINE = "HEK293T"
GROUP = "B_interoceptive"
SEED_SALT = "homeostatic-locality-v2"


def machinery(module):
    g = set()
    for s in ("core", "sensor", "tf", "feedback"):
        g |= set(pct.get_genes(module, s))
    try:
        g |= set(pct.get_loop(module)["feedback_core"])
    except KeyError:
        pass
    return g


def category(module):
    return pct.get_module_info(module)["category"]


def has_feedback_core(module):
    try:
        return len(pct.get_loop(module)["feedback_core"]) > 0
    except KeyError:
        return False


def _builder():
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location("pb17", Path(__file__).with_name("17_xatlas_pseudobulk.py"))
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_PB = _builder()


def resolve_genes(symbols, gpos):
    """Catalog symbols -> column indices, using the same alias table as the builder."""
    out = []
    for s in symbols:
        r = _PB.resolve(s, gpos)
        if r is not None:
            out.append(gpos[r])
    return sorted(set(out))


def rng_for(*parts):
    h = hashlib.sha256("|".join([SEED_SALT, *map(str, parts)]).encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "little"))


def load_pb(path):
    d = np.load(Path(path) / "pseudobulk.npz", allow_pickle=False)
    return {k: d[k] for k in d.files}


def locality_rows(pb, max_fold=MAX_FOLD, exclude_targets=True, off_module=False):
    """One row per (knockdown, module) with L = log2(own mean z^2 / median matched-null mean z^2)."""
    genes = pb["genes"].astype(str)
    umi = pb["ctrl_mean_umi"]
    univ = np.where(umi >= MIN_GENE_CTRL_UMI)[0]
    gpos = {genes[i]: i for i in univ}
    dec_edges = np.quantile(np.log10(umi[univ]), np.linspace(0, 1, N_DECILES + 1)[1:-1])
    decile = np.full(len(genes), -1)
    decile[univ] = np.searchsorted(dec_edges, np.log10(umi[univ]))
    by_dec = {d: univ[decile[univ] == d] for d in range(N_DECILES)}

    mods = pct.list_modules()
    act = {m: resolve_genes(pct.get_genes(m, "activity"), gpos) for m in mods}
    kd_row = {k: i for i, k in enumerate(pb["kd"].astype(str))}

    rows = []
    for m in mods:
        for kd in sorted(machinery(m)):
            if kd not in kd_row:
                continue
            if exclude_targets and kd in set(pct.get_genes(m, "activity")):
                continue
            i = kd_row[kd]
            eff = (pb["n_cells"][i] >= MIN_CELLS and np.isfinite(pb["fold"][i])
                   and pb["kd_ctrl_mean_umi"][i] >= MIN_KD_CTRL_UMI and pb["fold"][i] <= max_fold)
            if not eff:
                continue
            read_m = m
            if off_module:
                own = set().union(*(act[o] for o in mods if kd in machinery(o)))
                others = [o for o in mods if kd not in machinery(o) and len(act[o]) >= MIN_ACTIVITY_GENES
                          and not (set(act[o]) & own)]
                read_m = others[rng_for("off", m, kd).integers(len(others))]
            kd_cols = set(resolve_genes([kd], gpos))
            A = [j for j in act[read_m] if j not in kd_cols]
            if len(A) < MIN_ACTIVITY_GENES:
                continue
            z2 = pb["Z"][i].astype(np.float64) ** 2
            own_z2 = z2[A].mean()
            excl = set(A) | kd_cols
            rng = rng_for(m, kd, read_m)
            pools = {}
            for j in A:
                d = decile[j]
                if d not in pools:
                    pools[d] = np.array([g for g in by_dec[d] if g not in excl])
            draws = np.stack([rng.choice(pools[decile[j]], N_NULL) for j in A], 1)
            null = z2[draws].mean(1)
            rows.append({"module": m, "kd": kd, "read_module": read_m, "n_cells": int(pb["n_cells"][i]),
                         "fold": float(pb["fold"][i]), "n_activity": len(A), "own_mean_z2": float(own_z2),
                         "null_median_z2": float(np.median(null)),
                         "L": float(np.log2(max(own_z2, 1e-12) / max(np.median(null), 1e-12)))})
    return rows, act, umi


def module_table(rows, act, umi, min_kd=1):
    out = []
    for m in pct.list_modules():
        L = [r["L"] for r in rows if r["module"] == m]
        if len(L) < min_kd or len(act[m]) < MIN_ACTIVITY_GENES:
            continue
        out.append({"module": m, "category": category(m), "feedback_core": has_feedback_core(m),
                    "n_kd": len(L), "L_module": float(np.median(L)),
                    "baseline": float(np.log10(np.mean(umi[act[m]]) + 1e-3))})
    return out


def contrast(mt, group_fn, adjust=True):
    if not mt:
        return {"evaluable": False, "n_group": 0, "n_rest": 0}
    y = np.array([r["L_module"] for r in mt])
    x = np.array([r["baseline"] for r in mt])
    if adjust and len(mt) >= 3:
        b = np.polyfit(x, y, 1)
        y = y - np.polyval(b, x)
        slope = float(b[0])
    else:
        slope = None
    g = np.array([group_fn(r) for r in mt])
    n1, n0 = int(g.sum()), int((~g).sum())
    res = {"n_group": n1, "n_rest": n0, "baseline_slope": slope,
           "evaluable": n1 >= MIN_PER_GROUP and n0 >= MIN_PER_GROUP}
    if n1 and n0:
        res["median_group"] = float(np.median(y[g]))
        res["median_rest"] = float(np.median(y[~g]))
        res["p_one_sided"] = float(stats.mannwhitneyu(y[g], y[~g], alternative="greater").pvalue)
        res["auc"] = float(stats.mannwhitneyu(y[g], y[~g]).statistic / (n1 * n0))
    return res


def analyse_line(pb):
    rows, act, umi = locality_rows(pb)
    mt = module_table(rows, act, umi)
    is_intero = lambda r: r["category"] == GROUP  # noqa: E731
    out = {"primary": contrast(mt, is_intero)}
    rows_off, _, _ = locality_rows(pb, off_module=True)
    out["N1_off_module"] = contrast(module_table(rows_off, act, umi), is_intero)
    out["S1_feedback_core"] = contrast(mt, lambda r: r["feedback_core"])
    out["S2_unadjusted"] = contrast(mt, is_intero, adjust=False)
    rows_t, _, _ = locality_rows(pb, exclude_targets=False)
    out["S3_with_target_kds"] = contrast(module_table(rows_t, act, umi), is_intero)
    out["S4_min2_kd"] = contrast(module_table(rows, act, umi, min_kd=2), is_intero)
    rows_f, _, _ = locality_rows(pb, max_fold=SENS_MAX_FOLD)
    out["S5_fold_0.6"] = contrast(module_table(rows_f, act, umi), is_intero)
    return out, rows, mt


def verdict(res):
    p, r = res.get(PRIMARY_LINE), res.get(REPLICATION_LINE)

    def sig(x):
        return x is not None and x["primary"]["evaluable"] and x["primary"].get("p_one_sided", 1) < ALPHA

    def n1_violated(x):
        return x is not None and x["N1_off_module"].get("p_one_sided", 1) < ALPHA

    if p is None or not p["primary"]["evaluable"]:
        v = "NOT_EVALUABLE"
    elif not sig(p):
        v = "NULL"
    elif sig(r):
        v = "CONFIRMED"
    else:
        v = "PARTIAL"
    flags = [ln for ln in (PRIMARY_LINE, REPLICATION_LINE) if n1_violated(res.get(ln))]
    if flags and v in ("CONFIRMED", "PARTIAL"):
        v = v + "_NONSPECIFIC"
    return {"verdict": v, "N1_violated_in": flags,
            "primary_p": p["primary"].get("p_one_sided") if p else None,
            "replication_p": r["primary"].get("p_one_sided") if r else None,
            "replication_evaluable": bool(r and r["primary"]["evaluable"])}


def write_tsv(path, rows):
    if not rows:
        Path(path).write_text("")
        return
    cols = list(rows[0])
    with open(path, "w") as f:
        w = csv.DictWriter(f, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", action="append", required=True, help="NAME=dir with pseudobulk.npz")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    for spec in a.line:
        name, d = spec.split("=", 1)
        r, rows, mt = analyse_line(load_pb(d))
        res[name] = r
        write_tsv(out / f"{name}_per_kd.tsv", rows)
        write_tsv(out / f"{name}_per_module.tsv", mt)
    v = verdict(res)
    (out / "results.json").write_text(json.dumps(res, indent=1))
    (out / "verdict.json").write_text(json.dumps(v, indent=1))
    print(json.dumps(v, indent=1))


if __name__ == "__main__":
    main()
