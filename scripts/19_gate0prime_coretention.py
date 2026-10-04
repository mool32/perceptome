"""Gate 0' — are perceptome modules co-retained beyond the phylogenetic breadth of their own genes?

Pre-registered in docs/PREREG_gate0prime.md (sha256 in docs/PREREG_LOCK_gate0prime.txt).
Do not edit after the lock.

Why: in Paper 7's P5 the null (per-species coverage shuffle) is nearly the same for every module, so
the module z is ~ the number of species in which the module is intact, and module age is read from the
same presence row. "Old modules are coherent" is then close to a definition. Here each gene of a module
is replaced, in the null, by a random gene of the same age (deepest species with an ortholog) and
phylogenetic breadth (number of species with an ortholog), so gene histories are held fixed and only co-retention of the module's genes is tested.

Inputs (written by the operator from the Paper 7 pipeline, see the instruction):
  --species   TSV: species, divergence_mya        (the 15 species of P5)
  --presence  TSV: gene, <species columns> in {0,1} (gene has an ortholog, same rule as step04)
  --catalog   NAME=modules.json   ({module: [genes]}), repeatable; first one is the reference C0
  --check     NAME=per_module.tsv (pipeline output: module, age_mya, observed_intact), optional, repeatable
Usage:
  python scripts/19_gate0prime_coretention.py --species species.tsv --presence presence.tsv \
     --catalog C0=C0.json --catalog C2=C2.json ... --check C0=per_module_C0.tsv ... --out results/gate0prime
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats

# ---- locked parameters ---------------------------------------------------------
K = 0.6                  # module intact in a species if >= K of its genes have an ortholog there
N_NULL = 10000
MIN_POOL = 5             # widen the breadth window only if the exact (age, breadth) class has < 5 genes
ANCIENT_MYA = 1500
YOUNG_MYA = 800
MIN_PER_GROUP = 5
ALPHA = 0.05
SEED = 20261004
MAX_WIDENED_FRAC = 0.10  # catalog not evaluable if more of its genes needed a widened null class
CALIB_REPS = 3           # pseudo-modules per real module for the calibration check
MAX_CALIB_BIAS = 0.20    # |mean z_coh| of pseudo-modules above this -> catalog not evaluable
REFERENCE = "C0"


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_presence(path, species):
    rows = read_tsv(path)
    genes = np.array([r["gene"] for r in rows])
    P = np.array([[int(float(r[s])) for s in species] for r in rows], dtype=np.int8)
    return genes, P


def intact(Pm):
    """Pm: genes x species presence of one module -> boolean per species."""
    return Pm.mean(0) >= K - 1e-9


def module_age(intact_row, div):
    return float(div[intact_row].max()) if intact_row.any() else 0.0


def _pools(idx, breadth, gage):
    """Null pool per gene: same gene age (deepest species with an ortholog), breadth within +-w,
    w widened until >= MIN_POOL genes; if the age class is too small, the age constraint is dropped."""
    own = np.zeros(len(breadth), dtype=bool)
    own[idx] = True
    pools, widened = [], 0
    for i in idx:
        same = (gage == gage[i]) & ~own
        base = same if same.sum() >= MIN_POOL else ~own
        w = 0
        while True:
            cand = np.where(base & (np.abs(breadth - breadth[i]) <= w))[0]
            if len(cand) >= MIN_POOL or w > breadth.max():
                break
            w += 1
        widened += int(w > 0 or base is not same)
        pools.append(cand)
    return pools, widened


def _draw_distinct(pools, rng, n_null):
    """n_null rows, one gene per module gene from its pool, no gene twice in a row.
    Genes sharing an identical pool get distinct draws by a random ordering of that pool; rare
    collisions between different (widened) pools are redrawn."""
    draws = np.empty((n_null, len(pools)), dtype=np.int64)
    groups = {}
    for j, p in enumerate(pools):
        groups.setdefault(p.tobytes(), []).append(j)
    for key, js in groups.items():
        p = pools[js[0]]
        if len(js) <= len(p):
            order = np.argsort(rng.random((n_null, len(p))), 1)[:, :len(js)]
            draws[:, js] = p[order]
        else:                                   # pool smaller than the group: duplicates unavoidable
            draws[:, js] = rng.choice(p, (n_null, len(js)))
    for _ in range(100):
        srt = np.sort(draws, 1)
        bad = np.where((srt[:, 1:] == srt[:, :-1]).any(1))[0]
        if len(bad) == 0:
            break
        for r in bad:
            seen = set()
            for j, p in enumerate(pools):
                if draws[r, j] in seen:
                    free = np.setdiff1d(p, list(seen))
                    if len(free):
                        draws[r, j] = rng.choice(free)
                seen.add(int(draws[r, j]))
    return draws


def _midp_z(obs, null):
    q = (null < obs - 1e-12).mean() + 0.5 * (np.abs(null - obs) <= 1e-12).mean()
    q = min(max(q, 0.5 / len(null)), 1 - 0.5 / len(null))
    return float(q), float(stats.norm.ppf(q))


def coherence(genes_m, gidx, P, breadth, gage, rng, n_null=N_NULL):
    """Co-retention beyond gene breadth.

    Primary statistic S = variance over species of the module's coverage (fraction of its genes with an
    ortholog). With every gene's breadth fixed, S grows when the genes are gained and lost together and
    shrinks when they replace each other. Null: each gene replaced by a random gene of the same age
    and breadth.
    Secondary: number of species where the module is intact (>= K), against the same null.
    """
    idx = np.array([gidx[g] for g in genes_m])
    pools, widened = _pools(idx, breadth, gage)
    draws = _draw_distinct(pools, rng, n_null)                              # n_null x n_genes
    cov_null = P[draws].mean(1)                                              # n_null x species
    cov_obs = P[idx].mean(0)
    S_obs, S_null = float(cov_obs.var()), cov_null.var(1)
    q, z = _midp_z(S_obs, S_null)
    obs_int = int((cov_obs >= K - 1e-9).sum())
    int_null = (cov_null >= K - 1e-9).sum(1)
    _, z_int = _midp_z(obs_int, int_null)
    return {"observed_intact": obs_int, "S_obs": S_obs, "S_null_mean": float(S_null.mean()),
            "q": q, "z_coh": z, "z_intact_matched": z_int, "n_widened": widened}


def gene_ages(P, div):
    return np.array([div[r.astype(bool)].max() if r.any() else 0.0 for r in P])


def analyse_catalog(mods, genes, P, div, breadth, rng):
    gidx = {g: i for i, g in enumerate(genes)}
    gage = gene_ages(P, div)
    rows = []
    for m, gl in sorted(mods.items()):
        gl = [g for g in dict.fromkeys(gl) if g in gidx]
        if len(gl) < 2:
            continue
        idx = np.array([gidx[g] for g in gl])
        r = {"module": m, "n_genes": len(gl), "age_mya": module_age(intact(P[idx]), div)}
        r.update(coherence(gl, gidx, P, breadth, gage, rng))
        rows.append(r)
    return rows


def calibration(mods, genes, P, div, breadth, rng):
    """Pseudo-modules: each real module's genes replaced by random genes of the same gene age.
    Under no co-retention their z_coh must average ~0; returns the mean and SE."""
    gidx = {g: i for i, g in enumerate(genes)}
    gage = gene_ages(P, div)
    zs = []
    for _ in range(CALIB_REPS):
        for m, gl in sorted(mods.items()):
            gl = [g for g in dict.fromkeys(gl) if g in gidx]
            if len(gl) < 2:
                continue
            pick = set()
            for g in gl:
                cand = np.setdiff1d(np.where(gage == gage[gidx[g]])[0], list(pick))
                pick.add(int(rng.choice(cand)))
            zs.append(coherence(list(genes[sorted(pick)]), gidx, P, breadth, gage, rng)["z_coh"])
    zs = np.array(zs)
    return {"calib_mean_z": float(zs.mean()), "calib_se": float(zs.std() / np.sqrt(len(zs))), "calib_n": int(len(zs))}


def criteria(rows, calib=None):
    z = np.array([r["z_coh"] for r in rows])
    age = np.array([r["age_mya"] for r in rows])
    anc, yng = z[age >= ANCIENT_MYA], z[age <= YOUNG_MYA]
    out = {"n_ancient": int(len(anc)), "n_young": int(len(yng)),
           "median_z_ancient": float(np.median(anc)) if len(anc) else None,
           "median_z_young": float(np.median(yng)) if len(yng) else None}
    n_genes = sum(r["n_genes"] for r in rows)
    out["widened_frac"] = sum(r["n_widened"] for r in rows) / max(n_genes, 1)
    if calib:
        out.update(calib)
    bad_null = out["widened_frac"] > MAX_WIDENED_FRAC or (calib and abs(calib["calib_mean_z"]) > MAX_CALIB_BIAS)
    out["null_ok"] = not bad_null
    if len(anc) < MIN_PER_GROUP or len(yng) < MIN_PER_GROUP or bad_null:
        out.update(evaluable=False, L1=None, L2=None, L3=None, two_layer=None, ancient_units=None)
        return out
    p1 = float(stats.wilcoxon(anc, alternative="greater").pvalue) if np.any(anc != 0) else 1.0
    p2 = float(stats.wilcoxon(yng, alternative="less").pvalue) if np.any(yng != 0) else 1.0
    p3 = float(stats.mannwhitneyu(anc, yng, alternative="greater").pvalue)
    out.update(evaluable=True, p_L1_ancient_coherent=p1, p_L2_young_anticoherent=p2, p_L3_ancient_gt_young=p3,
               L1=p1 < ALPHA, L2=p2 < ALPHA, L3=p3 < ALPHA)
    out["two_layer"] = bool(out["L1"] and out["L2"] and out["L3"])
    out["ancient_units"] = bool(out["L1"] and out["L3"])
    rho = stats.spearmanr(age, z)
    out["spearman_age_zcoh"] = [float(rho.statistic), float(rho.pvalue)]
    return out


def verdict(per_cat, alternatives):
    ev = [c for c in alternatives if per_cat[c]["evaluable"]]
    if len(ev) < 2:
        return "NOT_EVALUABLE"
    if sum(per_cat[c]["two_layer"] for c in ev) >= 2:
        return "PASS_TWO_LAYERS"
    if sum(per_cat[c]["ancient_units"] for c in ev) >= 2:
        return "PASS_ANCIENT_UNITS_ONLY"
    return "FAIL"


def check(rows, path):
    """Recomputed intact counts and ages must equal the pipeline's for every module."""
    ref = {r["module"]: r for r in read_tsv(path)}
    bad = []
    for r in rows:
        p = ref.get(r["module"])
        if p is None:
            bad.append((r["module"], "missing in pipeline table"))
            continue
        if int(float(p["observed_intact"])) != r["observed_intact"] or abs(float(p["age_mya"]) - r["age_mya"]) > 0.5:
            bad.append((r["module"], f"pipeline {p['observed_intact']}/{p['age_mya']} vs {r['observed_intact']}/{r['age_mya']}"))
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--presence", required=True)
    ap.add_argument("--catalog", action="append", required=True)
    ap.add_argument("--check", action="append", default=[])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sp = read_tsv(a.species)
    species = [r["species"] for r in sp]
    div = np.array([float(r["divergence_mya"]) for r in sp])
    genes, P = load_presence(a.presence, species)
    breadth = P.sum(1)
    checks = dict(c.split("=", 1) for c in a.check)
    per_cat, mismatches = {}, {}
    names = []
    for spec in a.catalog:
        name, path = spec.split("=", 1)
        names.append(name)
        mods = json.load(open(path))
        rng = np.random.default_rng([SEED, sum(map(ord, name))])
        rows = analyse_catalog(mods, genes, P, div, breadth, rng)
        calib = calibration(mods, genes, P, div, breadth, rng)
        with open(out / f"zcoh_{name}.tsv", "w") as f:
            w = csv.DictWriter(f, list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        if name in checks:
            mismatches[name] = check(rows, checks[name])
        per_cat[name] = criteria(rows, calib)
    alternatives = [n for n in names if n != REFERENCE]
    res = {"verdict": verdict(per_cat, alternatives), "reference": REFERENCE, "alternatives": alternatives,
           "per_catalog": per_cat, "check_mismatches": mismatches}
    if any(mismatches.values()):
        res["verdict"] = "STOP_INPUT_MISMATCH"
    (out / "gate0prime.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({"verdict": res["verdict"], **{k: {kk: v.get(kk) for kk in ("two_layer", "ancient_units", "n_ancient", "n_young")}
                                                     for k, v in per_cat.items()}}, indent=1))


if __name__ == "__main__":
    main()
