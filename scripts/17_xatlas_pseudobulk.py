"""Pseudobulk builder for the homeostatic-locality pre-registration v2 (X-Atlas/Orion).

Locked together with docs/PREREG_homeostatic_locality_v2.md (sha256 in
docs/PREREG_LOCK_v2.txt). Do not edit after the lock; deviations go to
docs/PREREG_homeostatic_locality_v2_DEVIATIONS.md.

Reads raw UMI counts only for (i) a seeded sample of control cells and (ii) cells whose
target is one of the knockdowns the analysis can use (module machinery genes, see
`analysis_knockdowns`). Everything else in the 8M-cell atlas is never read.

Per knockdown it writes the mean per-cell z-score of log-CP10k against control cells of the
same batch (Replogle-style), the knockdown fold of the target gene, and the cell count.

Two modes:
  --inspect   prints obs/var column names, dtypes, X layout and the labels that look like
              controls. Reads no expression value. Run this first to set the column flags.
  (default)   builds the pseudobulk.

The input may be a local path or an fsspec URL (e.g. hf://datasets/<org>/<repo>/<file>),
read lazily with range requests, so the full file never has to be on disk.

Usage:
  python scripts/17_xatlas_pseudobulk.py --h5ad <file-or-url> --inspect
  python scripts/17_xatlas_pseudobulk.py --h5ad <file-or-url> --out <dir> \
      --target_col gene_target --control_label <label> [--batch_col <col>] \
      [--gene_col <var column with symbols>] [--filter_col pass_guide_filter]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perceptome as pct  # noqa: E402

# ---- locked parameters --------------------------------------------------------
SEED = 0
N_CTRL = 10000             # control cells sampled for mu / sigma (all if fewer)
MIN_CTRL_PER_BATCH = 200   # fewer -> that batch uses the pooled control mu / sigma
SIGMA_FLOOR = 0.05         # floor on control SD of log-CP10k (avoids the +inf of Replogle's file)
CHUNK_ROWS = 4000

# Symbol aliases: catalog symbol -> other names the atlas may use (HGNC renames).
ALIASES = {
    "TMEM173": ["STING1"], "STING1": ["TMEM173"],
    "ARNTL": ["BMAL1"], "BMAL1": ["ARNTL"],
    "MB21D1": ["CGAS"], "CGAS": ["MB21D1"],
    "G6PC": ["G6PC1"], "ACPP": ["ACP3"], "CTGF": ["CCN2"], "CYR61": ["CCN1"],
    "PHD2": ["EGLN1"], "GILZ": ["TSC22D3"], "COP1": ["RFWD2"], "RFWD2": ["COP1"],
    "KIAA0101": ["PCLAF"], "FAM129A": ["NIBAN1"], "C10orf10": ["DEPP1"],
}


def machinery(module):
    """Genes whose knockdown counts as perturbing the module (same rule as the exploration)."""
    g = set()
    for s in ("core", "sensor", "tf", "feedback"):
        g |= set(pct.get_genes(module, s))
    try:
        g |= set(pct.get_loop(module)["feedback_core"])
    except KeyError:
        pass
    return g


def analysis_knockdowns():
    """Catalog symbols of every knockdown the locked analysis can use."""
    out = set()
    for m in pct.list_modules():
        out |= machinery(m)
    return sorted(out)


def catalog_genes():
    out = set(analysis_knockdowns())
    for m in pct.list_modules():
        out |= set(pct.get_genes(m, "activity"))
    return out


def resolve(symbol, available):
    if symbol in available:
        return symbol
    for a in ALIASES.get(symbol, []):
        if a in available:
            return a
    return None


# ---- h5ad access ---------------------------------------------------------------
def open_h5(path):
    import h5py
    if "://" in path:
        import fsspec
        f = fsspec.open(path, "rb", block_size=16 * 2 ** 20, cache_type="blockcache").open()
        return h5py.File(f, "r")
    return h5py.File(path, "r")


def read_column(group, name):
    """Read an obs/var column, decoding anndata categoricals and byte strings."""
    import h5py
    node = group[name]
    if isinstance(node, h5py.Group) and "categories" in node:
        cats = np.asarray(node["categories"][()]).astype(str)
        codes = np.asarray(node["codes"][()])
        out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "")
        return out.astype(str)
    v = np.asarray(node[()])
    return v.astype(str) if v.dtype.kind in "SO" else v


def index_column(group):
    name = group.attrs.get("_index", "_index")
    return read_column(group, name)


class Rows:
    """Row reader for X stored dense or as CSR; returns float32 dense blocks for gene columns `cols`."""

    def __init__(self, h, cols):
        import h5py
        X = h["X"]
        self.cols = np.asarray(cols)
        self.dense = isinstance(X, h5py.Dataset)
        self.X = X
        if not self.dense:
            enc = X.attrs.get("encoding-type", "csr_matrix")
            if "csr" not in str(enc):
                raise SystemExit(f"X is {enc}; only CSR or dense is supported")
            self.indptr = np.asarray(X["indptr"][()], dtype=np.int64)
            self.colmap = np.full(int(X.attrs["shape"][1]), -1, dtype=np.int64)
            self.colmap[self.cols] = np.arange(len(self.cols))

    def read(self, idx):
        idx = np.sort(np.asarray(idx, dtype=np.int64))
        out = np.zeros((len(idx), len(self.cols)), dtype=np.float32)
        tot = np.zeros(len(idx))
        if self.dense:
            for s in range(0, len(idx), CHUNK_ROWS // 4):
                b = idx[s:s + CHUNK_ROWS // 4]
                blk = np.asarray(self.X[b])
                tot[s:s + len(b)] = blk.sum(1)
                out[s:s + len(b)] = blk[:, self.cols]
            return idx, out, tot
        # coalesce runs of nearby rows into one contiguous read of data/indices
        start = 0
        while start < len(idx):
            end = start + 1
            while (end < len(idx) and idx[end] - idx[end - 1] <= 64 and end - start < CHUNK_ROWS):
                end += 1
            r0, r1 = idx[start], idx[end - 1] + 1
            p0, p1 = self.indptr[r0], self.indptr[r1]
            data = np.asarray(self.X["data"][p0:p1], dtype=np.float32)
            ind = np.asarray(self.X["indices"][p0:p1], dtype=np.int64)
            for k in range(start, end):
                r = idx[k]
                a, b = self.indptr[r] - p0, self.indptr[r + 1] - p0
                tot[k] = data[a:b].sum()
                c = self.colmap[ind[a:b]]
                keep = c >= 0
                out[k, c[keep]] = data[a:b][keep]
            start = end
        return idx, out, tot


def lognorm(counts, totals):
    return np.log1p(counts / np.maximum(totals[:, None], 1.0) * 1e4)


# ---- modes ---------------------------------------------------------------------
def inspect(path):
    import h5py
    h = open_h5(path)
    obs, var = h["obs"], h["var"]
    info = {"obs_columns": sorted(k for k in obs.keys() if k != "__categories"),
            "var_columns": sorted(k for k in var.keys() if k != "__categories"),
            "X_layout": "dense" if isinstance(h["X"], h5py.Dataset) else dict(h["X"].attrs).get("encoding-type"),
            "X_shape": list(h["X"].shape) if isinstance(h["X"], h5py.Dataset) else [int(x) for x in h["X"].attrs["shape"]],
            "X_dtype": str(h["X"].dtype if isinstance(h["X"], h5py.Dataset) else h["X"]["data"].dtype),
            "layers": sorted(h["layers"].keys()) if "layers" in h else []}
    for col in info["obs_columns"]:
        node = obs[col]
        if isinstance(node, h5py.Group) and "categories" in node:
            cats = np.asarray(node["categories"][()]).astype(str)
            info.setdefault("obs_categorical_n", {})[col] = len(cats)
            hits = [c for c in cats if any(t in c.lower() for t in ("non", "control", "ntc", "safe", "scram"))]
            if hits:
                info.setdefault("control_like_labels", {})[col] = hits[:20]
    for col in info["var_columns"]:
        v = read_column(var, col)
        info.setdefault("var_examples", {})[col] = [str(x) for x in v[:5]]
    info["var_index_examples"] = [str(x) for x in index_column(var)[:5]]
    print(json.dumps(info, indent=1))


def build(a):
    t0 = time.time()
    h = open_h5(a.h5ad)
    obs, var = h["obs"], h["var"]
    target = read_column(obs, a.target_col).astype(str)
    keep = np.ones(len(target), dtype=bool)
    if a.filter_col:
        keep &= read_column(obs, a.filter_col).astype(bool)
    batch = read_column(obs, a.batch_col).astype(str) if a.batch_col else np.full(len(target), "all")
    genes = (read_column(var, a.gene_col) if a.gene_col else index_column(var)).astype(str)
    gset = set(genes)
    n_genes = len(genes)

    # readout universe: every measured gene (the analysis restricts it further)
    cols = np.arange(n_genes)
    rows = Rows(h, cols)

    rng = np.random.default_rng(SEED)
    ctrl_all = np.where((target == a.control_label) & keep)[0]
    if len(ctrl_all) == 0:
        raise SystemExit(f"no control cells with label {a.control_label!r}")
    ctrl = np.sort(rng.choice(ctrl_all, min(N_CTRL, len(ctrl_all)), replace=False))

    labels_avail = set(np.unique(target[keep]))
    wanted = {}
    for s in analysis_knockdowns():
        r = resolve(s, labels_avail)
        if r is not None:
            wanted[s] = r

    print(f"[prep] cells {len(target)}, controls {len(ctrl_all)} (sampled {len(ctrl)}), "
          f"knockdowns wanted {len(analysis_knockdowns())}, present {len(wanted)}", flush=True)

    # controls: per-batch mu / sigma of log-CP10k, mean UMI per gene
    cidx, C, ctot = rows.read(ctrl)
    ctrl_mean_umi = C.mean(0)
    ctrl_mean_cp10k = (C / np.maximum(ctot[:, None], 1.0) * 1e4).mean(0)
    L = lognorm(C, ctot)
    del C
    cb = batch[cidx]
    pooled_mu, pooled_sd = L.mean(0), np.maximum(L.std(0), SIGMA_FLOOR)
    stats_b = {}
    for b in np.unique(cb):
        m = cb == b
        if m.sum() >= MIN_CTRL_PER_BATCH:
            stats_b[b] = (L[m].mean(0), np.maximum(L[m].std(0), SIGMA_FLOOR))
    ctrl_mean_log = pooled_mu
    del L
    print(f"[prep] controls done ({time.time() - t0:.0f}s); batches with own stats: {len(stats_b)}", flush=True)

    kd_names, Z, ncell, fold, kd_ctrl_umi = [], [], [], [], []
    for i, (sym, lab) in enumerate(sorted(wanted.items())):
        cells = np.where((target == lab) & keep)[0]
        if len(cells) == 0:
            continue
        idx, K, tot = rows.read(cells)
        Lk = lognorm(K, tot)
        zsum = np.zeros(n_genes)
        for b in np.unique(batch[idx]):
            m = batch[idx] == b
            mu, sd = stats_b.get(b, (pooled_mu, pooled_sd))
            zsum += ((Lk[m] - mu) / sd).sum(0)
        g = resolve(sym, gset)
        if g is not None:
            gi = int(np.where(genes == g)[0][0])
            cp_kd = (K[:, gi] / np.maximum(tot, 1) * 1e4).mean()
            cp_c = ctrl_mean_cp10k[gi]
            f = float(cp_kd / cp_c) if cp_c > 0 else float("nan")
            cu = float(ctrl_mean_umi[gi])
        else:
            f, cu = float("nan"), float("nan")
        kd_names.append(sym)
        Z.append((zsum / len(idx)).astype(np.float32))
        ncell.append(len(idx))
        fold.append(f)
        kd_ctrl_umi.append(cu)
        if (i + 1) % 25 == 0:
            print(f"[prep] {i + 1}/{len(wanted)} knockdowns ({time.time() - t0:.0f}s)", flush=True)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "pseudobulk.npz", Z=np.vstack(Z), genes=genes, kd=np.array(kd_names),
                        n_cells=np.array(ncell), fold=np.array(fold), kd_ctrl_mean_umi=np.array(kd_ctrl_umi),
                        ctrl_mean_umi=ctrl_mean_umi, ctrl_mean_log=ctrl_mean_log)
    meta = {"h5ad": a.h5ad, "target_col": a.target_col, "control_label": a.control_label,
            "batch_col": a.batch_col, "gene_col": a.gene_col, "filter_col": a.filter_col,
            "n_cells_total": int(len(target)), "n_cells_kept": int(keep.sum()),
            "n_controls": int(len(ctrl_all)), "n_controls_sampled": int(len(ctrl)),
            "n_batches_own_stats": len(stats_b), "n_genes": int(n_genes),
            "n_knockdowns_wanted": len(analysis_knockdowns()), "n_knockdowns_written": len(kd_names),
            "missing_knockdowns": sorted(set(analysis_knockdowns()) - set(wanted)),
            "catalog_genes_missing_from_var": sorted(g for g in catalog_genes() if resolve(g, gset) is None),
            "runtime_s": round(time.time() - t0, 1)}
    (out / "pseudobulk_meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps({k: v for k, v in meta.items() if not isinstance(v, list)}, indent=1))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--inspect", action="store_true")
    p.add_argument("--out")
    p.add_argument("--target_col", default="gene_target")
    p.add_argument("--control_label")
    p.add_argument("--batch_col")
    p.add_argument("--gene_col")
    p.add_argument("--filter_col")
    a = p.parse_args()
    if a.inspect:
        inspect(a.h5ad)
    else:
        if not (a.out and a.control_label):
            p.error("--out and --control_label are required")
        build(a)


if __name__ == "__main__":
    main()
