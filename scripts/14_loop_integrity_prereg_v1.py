"""Loop-integrity test, pre-registration v1 (docs/PREREG_loop_integrity_v1.md).

This script IS the pre-registered analysis: it is locked together with the
pre-registration (sha256 in docs/PREREG_LOCK.txt). Do not edit it after the lock;
any change is a deviation and goes to docs/PREREG_loop_integrity_v1_DEVIATIONS.md.

Usage:
  python scripts/14_loop_integrity_prereg_v1.py \
      --h5ad K562_gwps_raw_singlecell_01.h5ad \
      --out results/loop_integrity_v1

Inputs: raw UMI counts (Replogle et al. 2022, K562 genome-wide CRISPRi), the
perturbation table data/evidence/replogle2022_K562_gwps_contrasts.tsv, and the
loop layer (pct.loop_output_genes / pct.load_loops).
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perceptome as pct  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CONTRASTS = ROOT / "data" / "evidence" / "replogle2022_K562_gwps_contrasts.tsv"

# ---- locked parameters (§3-§7) ------------------------------------------------
SEED = 0
MIN_CELLS = 50
MAX_FOLD = 0.50                 # author fold_expr <= 0.50  (>= 50 % knockdown)
SENS_MAX_FOLD = 0.60            # UPR-ATF6 sensitivity analysis only
CTRL_PER_CELL = 5
DEPTH_PCT = 10                  # D = 10th percentile of total UMI over selected controls
MEAN_MIN = 0.1                  # gene enters M2 if mean >= 0.1 in both groups
PHI_EPS = 0.01
LFC_EPS = 0.1
MIN_GENES = 3                   # a (KD, module) pair needs >= 3 output genes
N_NULL = 1000
N_BOOT = 500
N_PERM = 10000
N_PSEUDO = 50
PSEUDO_CELLS = 150
ALPHA = 0.05
Z_CRIT = 1.96
N1_Z_CRIT = 2.128               # one-sided, Bonferroni over the 3 negative-control break KDs (0.05/3)
G2_MAX_FRAC = 0.10
MIN_LOOPS_G1 = 3

PRIMARY = ["ERK/MAPK", "HIF", "JAK-STAT", "SREBP", "UPR-IRE1", "UPR-PERK"]
SECONDARY_BREAK_ONLY = ["NRF2", "HSF1"]
INPUT_UP_SETS = {"mTOR": ["TSC1"], "ERK/MAPK": ["NF1"]}
NEGATIVE_CONTROL = "p53"
SENSITIVITY = "UPR-ATF6"
CONTROL_LABEL = "non-targeting"


# ---- perturbation table -------------------------------------------------------
def read_contrasts(path=CONTRASTS):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["n_cells"] = int(r["n_cells"])
            r["fold_expr"] = float(r["fold_expr"]) if r["fold_expr"] not in ("", None) else None
            rows.append(r)
    return rows


def eligible(r, max_fold=MAX_FOLD):
    return (r["in_screen"] == "yes" and r["n_cells"] >= MIN_CELLS
            and r["fold_expr"] is not None and r["fold_expr"] <= max_fold)


def analysis_set(rows):
    """{module: {class: [genes]}} for every module, eligibility applied."""
    out = {}
    for r in rows:
        mf = SENS_MAX_FOLD if r["module"] == SENSITIVITY else MAX_FOLD
        if eligible(r, mf):
            out.setdefault(r["module"], {}).setdefault(r["class"], []).append(r["gene"])
    return out


# ---- data access -------------------------------------------------------------
def open_screen(h5ad):
    import anndata as ad
    a = ad.read_h5ad(h5ad, backed="r")
    obs = a.obs
    gcol = next((c for c in ("gene", "perturbation") if c in obs), None)
    bcol = next((c for c in ("gem_group", "batch") if c in obs), None)
    if gcol is None or bcol is None:
        raise SystemExit(f"obs lacks gene/gem-group columns: {list(obs.columns)}")
    genes = (a.var["gene_name"].astype(str).to_numpy() if "gene_name" in a.var
             else np.asarray(a.var_names, dtype=str))
    return a, obs[gcol].astype(str).to_numpy(), obs[bcol].astype(str).to_numpy(), genes


def _rows(a, idx):
    X = a.X[np.asarray(idx)]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return np.rint(X).astype(np.int64)


def load_thinned(a, idx, ctrl_idx, rng, chunk=4000):
    """Read the needed rows in chunks and thin them to depth D on the fly.

    D = DEPTH_PCT-th percentile of total UMI over `ctrl_idx` (first pass reads
    control rows only for their totals). Returns (kept cell indices, thinned
    int16 counts, D, n dropped below D).
    """
    ctrl_idx = np.asarray(sorted(set(int(i) for i in ctrl_idx)))
    tot = np.concatenate([_rows(a, ctrl_idx[i:i + chunk]).sum(1) for i in range(0, len(ctrl_idx), chunk)])
    depth = int(np.percentile(tot, DEPTH_PCT))
    idx = np.asarray(sorted(set(int(i) for i in idx)))
    kept, mats = [], []
    for i in range(0, len(idx), chunk):
        part = idx[i:i + chunk]
        Xt, keep = thin(_rows(a, part), depth, rng)
        kept.append(part[keep])
        mats.append(Xt.astype(np.int16))
    return np.concatenate(kept), np.vstack(mats), depth, int(len(idx) - sum(len(k) for k in kept))


# ---- statistics ---------------------------------------------------------------
def thin(X, depth, rng):
    """Binomially thin each row to `depth` total counts; rows below depth are dropped."""
    tot = X.sum(1)
    keep = tot >= depth
    Xk, tk = X[keep], tot[keep]
    p = depth / tk
    return rng.binomial(Xk, p[:, None]), keep


def gene_stats(X):
    mu = X.mean(0)
    var = X.var(0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(mu > 0, (var - mu) / mu ** 2, np.nan)
    return mu, phi


def m1(mu_k, mu_c, gi):
    return float(np.mean(np.log2((mu_k[gi] + LFC_EPS) / (mu_c[gi] + LFC_EPS))))


def dvec(phi_k, phi_c):
    a = np.clip(phi_k, -PHI_EPS / 2, None) + PHI_EPS
    b = np.clip(phi_c, -PHI_EPS / 2, None) + PHI_EPS
    return np.log(a / b)


def compare(P, C, out_idx, pool_idx, rng, feedback_idx=None, boot=True):
    """All metrics for one perturbed group P vs controls C on one output gene set."""
    P = P.astype(np.float64)
    C = C.astype(np.float64)
    mu_k, phi_k = gene_stats(P)
    mu_c, phi_c = gene_stats(C)
    res = {"n_pert": int(P.shape[0]), "n_ctrl": int(C.shape[0]), "n_out": int(len(out_idx))}
    if len(out_idx) < MIN_GENES:
        res["status"] = "too_few_output_genes"
        return res
    res["M1"] = m1(mu_k, mu_c, out_idx)
    ok = (mu_k >= MEAN_MIN) & (mu_c >= MEAN_MIN) & np.isfinite(phi_k) & np.isfinite(phi_c)
    gi = np.array([g for g in out_idx if ok[g]], dtype=int)
    res["n_out_m2"] = int(len(gi))
    if len(gi) < MIN_GENES:
        res["status"] = "too_few_genes_for_M2"
        return res
    d = dvec(phi_k, phi_c)
    res["M2"] = float(np.median(d[gi]))
    # matched null: same size, same control-mean decile, pool genes passing the same filter
    pool = np.array([g for g in pool_idx if ok[g]], dtype=int)
    edges = np.quantile(mu_c[pool], np.linspace(0, 1, 11))
    dec = lambda x: np.clip(np.searchsorted(edges, x, side="right") - 1, 0, 9)  # noqa: E731
    pool_dec = dec(mu_c[pool])
    by_dec = {k: pool[pool_dec == k] for k in range(10)}
    need = dec(mu_c[gi])
    null = np.empty(N_NULL)
    for i in range(N_NULL):
        pick = [rng.choice(by_dec[k]) if len(by_dec[k]) else rng.choice(pool) for k in need]
        null[i] = np.median(d[np.array(pick)])
    sd = null.std(ddof=1)
    res["Z2"] = float((res["M2"] - null.mean()) / sd) if sd > 0 else float("nan")
    if boot:
        Ps, Cs = P[:, out_idx], C[:, out_idx]
        sub = np.arange(len(out_idx))
        gpos_in_sub = np.array([np.where(out_idx == g)[0][0] for g in gi])
        b1, b2 = np.empty(N_BOOT), np.empty(N_BOOT)
        for i in range(N_BOOT):
            Pb = Ps[rng.integers(0, Ps.shape[0], Ps.shape[0])]
            Cb = Cs[rng.integers(0, Cs.shape[0], Cs.shape[0])]
            mk, pk = gene_stats(Pb)
            mc, pc = gene_stats(Cb)
            b1[i] = m1(mk, mc, sub)
            b2[i] = np.nanmedian(dvec(pk, pc)[gpos_in_sub])
        res["M1_ci"] = [float(np.percentile(b1, 2.5)), float(np.percentile(b1, 97.5))]
        res["M1_ci99"] = [float(np.percentile(b1, 0.5)), float(np.percentile(b1, 99.5))]
        res["M2_ci"] = [float(np.nanpercentile(b2, 2.5)), float(np.nanpercentile(b2, 97.5))]
    if feedback_idx is not None and len(feedback_idx):
        def rho(X):
            L = np.log1p(X)
            z = (L - L.mean(0)) / np.where(L.std(0) > 0, L.std(0), 1)
            return stats.spearmanr(z[:, out_idx].mean(1), z[:, feedback_idx].mean(1)).statistic
        res["M3_delta_rho"] = float(rho(P) - rho(C))
    res["status"] = "ok"
    return res


def sign_flip_p(x, rng, n=N_PERM):
    x = np.asarray(x, float)
    obs = x.mean()
    flips = rng.choice([-1.0, 1.0], size=(n, len(x)))
    return float((np.sum((flips * np.abs(x)).mean(1) >= obs) + 1) / (n + 1))


def ols_class_coef(z, is_break, abs_m1):
    X = np.column_stack([np.ones(len(z)), is_break, abs_m1])
    beta, *_ = np.linalg.lstsq(X, z, rcond=None)
    resid = z - X @ beta
    dof = len(z) - X.shape[1]
    s2 = resid @ resid / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    t = beta[1] / se[1]
    return float(beta[1]), float(stats.t.sf(t, dof))


# ---- main ---------------------------------------------------------------------
def run(h5ad, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    a, glab, gem, var_genes = open_screen(h5ad)
    gpos = {g: i for i, g in enumerate(var_genes)}
    loops = pct.load_loops()["loops"]
    loop_genes = {g for e in loops.values() for k in ("feedback_core", "loop_targets") for g in e[k]}
    pool_idx = np.array([i for g, i in gpos.items() if g not in loop_genes], dtype=int)

    rows = read_contrasts()
    aset = analysis_set(rows)

    # (KD, module, class, set) units, frozen order
    units = []
    for m in PRIMARY:
        for cls in ("break", "input_down"):
            units += [(g, m, cls, "primary") for g in aset.get(m, {}).get(cls, [])]
    for m in SECONDARY_BREAK_ONLY:
        units += [(g, m, "break", "secondary_break_only") for g in aset.get(m, {}).get("break", [])]
    for m, gl in INPUT_UP_SETS.items():
        units += [(g, m, "input_up", "input_up") for g in gl if g in aset.get(m, {}).get("input_up", [])]
    units += [(g, "mTOR", "input_down", "input_up") for g in aset.get("mTOR", {}).get("input_down", [])]
    for cls in ("break", "input_down"):
        units += [(g, NEGATIVE_CONTROL, cls, "negative_control")
                  for g in aset.get(NEGATIVE_CONTROL, {}).get(cls, [])]
    for cls in ("break", "input_down"):
        units += [(g, SENSITIVITY, cls, "sensitivity") for g in aset.get(SENSITIVITY, {}).get(cls, [])]
    spec = [("MBTPS1", "SREBP"), ("MBTPS1", "UPR-ATF6"), ("MBTPS2", "SREBP"), ("MBTPS2", "UPR-ATF6"),
            ("ATF6", "SREBP"), ("ATF6", "UPR-ATF6"), ("SCAP", "SREBP"), ("SCAP", "UPR-ATF6")]
    units += [(g, m, "specificity", "specificity") for g, m in spec]

    # cell selection
    ctrl_all = np.where(glab == CONTROL_LABEL)[0]
    kd_cells = {g: np.where(glab == g)[0] for g in sorted({u[0] for u in units})}
    ctrl_sel = {}
    for g, cells in kd_cells.items():
        r = np.random.default_rng([SEED, sum(map(ord, g))])
        groups = set(gem[cells])
        pool = ctrl_all[np.isin(gem[ctrl_all], list(groups))]
        n = min(len(pool), CTRL_PER_CELL * len(cells))
        ctrl_sel[g] = np.sort(r.choice(pool, n, replace=False))
    pseudo = []
    pr = np.random.default_rng([SEED, 7])
    for _ in range(N_PSEUDO):
        p = np.sort(pr.choice(ctrl_all, PSEUDO_CELLS, replace=False))
        pool = np.setdiff1d(ctrl_all[np.isin(gem[ctrl_all], list(set(gem[p])))], p)
        c = np.sort(pr.choice(pool, min(len(pool), CTRL_PER_CELL * PSEUDO_CELLS), replace=False))
        pseudo.append((p, c))

    need = set(np.concatenate([*kd_cells.values(), *ctrl_sel.values(),
                               *[np.concatenate(x) for x in pseudo]]).tolist())
    ctrl_union = sorted(set(np.concatenate(list(ctrl_sel.values())).tolist()))
    kept_idx, Xt, depth, n_dropped = load_thinned(a, sorted(need), ctrl_union, rng)
    kept_row = {int(i): k for k, i in enumerate(kept_idx)}

    def mat(cells):
        r = [kept_row[int(i)] for i in cells if int(i) in kept_row]
        return Xt[r], len(cells) - len(r)

    def out_idx(module, kd):
        gs = [g for g in pct.loop_output_genes(module) if g != kd and g in gpos]
        return np.array([gpos[g] for g in gs], dtype=int), gs

    def fb_idx(module, kd):
        e = loops[module]
        if not e["feedback_readable"] or e["feedback_is_output"]:
            return None
        return np.array([gpos[g] for g in e["feedback_core"] if g != kd and g in gpos], dtype=int)

    results = []
    for kd, m, cls, setname in units:
        P, dropP = mat(kd_cells[kd])
        C, dropC = mat(ctrl_sel[kd])
        oi, og = out_idx(m, kd)
        r = compare(P, C, oi, pool_idx, rng, fb_idx(m, kd))
        r.update(kd=kd, module=m, cls=cls, set=setname, dropped_below_depth=[dropP, dropC], output_genes=og)
        results.append(r)

    # G2 pseudo-perturbations
    g2 = []
    for p, c in pseudo:
        P, _ = mat(p)
        C, _ = mat(c)
        for m in PRIMARY:
            oi, _ = out_idx(m, None)
            r = compare(P, C, oi, pool_idx, rng, boot=False)
            if "Z2" in r and np.isfinite(r["Z2"]):
                g2.append(abs(r["Z2"]) > Z_CRIT)
    g2_frac = float(np.mean(g2)) if g2 else float("nan")
    G2 = bool(g2 and g2_frac <= G2_MAX_FRAC)

    # G1 per primary loop
    prim = [r for r in results if r["set"] == "primary" and r.get("status") == "ok"]
    G1 = {}
    for m in PRIMARY:
        dn = [r for r in prim if r["module"] == m and r["cls"] == "input_down"]
        med = float(np.median([r["M1"] for r in dn])) if dn else float("nan")
        anyci = any(r["M1_ci"][1] < 0 for r in dn)
        G1[m] = {"median_M1_input_down": med, "any_ci_below_0": anyci, "pass": bool(dn and med < 0 and anyci)}
    loops_ok = [m for m in PRIMARY if G1[m]["pass"]]

    verdict = {"depth_D": depth, "n_cells_dropped_below_D": n_dropped, "G2": {"frac_abs_Z2_gt_1.96": g2_frac, "n": len(g2), "pass": G2},
               "G1": G1, "loops_passing_G1": loops_ok}
    if not G2:
        verdict["decision"] = "STOP: G2 failed (dispersion metric not distinguishable from sampling noise)"
    elif len(loops_ok) < MIN_LOOPS_G1:
        verdict["decision"] = f"STOP: fewer than {MIN_LOOPS_G1} loops pass G1 (readout too weak)"
    else:
        h = [r for r in prim if r["module"] in loops_ok and np.isfinite(r.get("Z2", np.nan))]
        zb = np.array([r["Z2"] for r in h if r["cls"] == "break"])
        zd = np.array([r["Z2"] for r in h if r["cls"] == "input_down"])
        p_flip = sign_flip_p(zb, np.random.default_rng([SEED, 1])) if len(zb) else 1.0
        p_mw = float(stats.mannwhitneyu(zb, zd, alternative="greater").pvalue) if len(zb) and len(zd) else 1.0
        H1 = bool(zb.mean() > 0 and p_flip < ALPHA and p_mw < ALPHA) if len(zb) else False
        coef, p_spec = ols_class_coef(np.array([r["Z2"] for r in h]),
                                      np.array([r["cls"] == "break" for r in h], float),
                                      np.array([abs(r["M1"]) for r in h]))
        H1spec = bool(coef > 0 and p_spec < ALPHA)
        nc = [r for r in results if r["set"] == "negative_control" and r.get("status") == "ok"]
        N1_violated = [r["kd"] for r in nc if r["cls"] == "break" and r.get("Z2", 0) >= N1_Z_CRIT]
        verdict.update({
            "H1": {"n_break": len(zb), "n_input_down": len(zd), "mean_Z2_break": float(zb.mean()) if len(zb) else None,
                   "p_sign_flip": p_flip, "p_mann_whitney": p_mw, "pass": H1},
            "H1_spec": {"break_coef": coef, "p_one_sided": p_spec, "pass": H1spec},
            "N1": {"violated_by": N1_violated, "not_violated": not N1_violated},
        })
        verdict["decision"] = ("POSITIVE: G2, H1, H1-spec pass and N1 not violated"
                               if (H1 and H1spec and not N1_violated) else "NULL: decision rule not met")

    (out / "units.json").write_text(json.dumps(results, indent=1))
    (out / "verdict.json").write_text(json.dumps(verdict, indent=1))
    with open(out / "units.tsv", "w") as f:
        cols = ["set", "module", "cls", "kd", "n_pert", "n_ctrl", "n_out", "n_out_m2",
                "M1", "M2", "Z2", "M3_delta_rho", "status"]
        f.write("\t".join(cols) + "\n")
        for r in results:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    return verdict


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    v = run(a.h5ad, a.out)
    print(json.dumps(v, indent=1))


if __name__ == "__main__":
    main()
