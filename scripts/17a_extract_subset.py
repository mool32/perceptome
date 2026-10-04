"""Sequential subset extractor for scripts/17_xatlas_pseudobulk.py (amendment 1 to pre-registration v2).

Why: script 17 reads the rows it needs by random access. Over HTTP each of the ~76k scattered rows
costs a whole cache block, so reading 10k controls alone moves ~160 GB and the connection drops.
This script makes ONE sequential pass over the remote CSR matrix (large range reads, retries,
resumable checkpoints) and keeps only the rows script 17 would read:

  - the control sample, chosen with exactly script 17's rule (same SEED, N_CTRL, filter), and
  - every cell of the knockdowns script 17 would use (same symbol resolution and filter).

It writes a small local h5ad with those rows, all genes, and the obs columns script 17 needs.
Script 17 then runs on it unchanged: with exactly N_CTRL controls present it samples all of them,
i.e. the same control set it would have drawn from the full file. No statistic is computed here.

Usage:
  python scripts/17a_extract_subset.py --h5ad <url-or-path> --out <subset.h5ad> --work <dir> \
      --target_col gene_target --control_label Non-Targeting \
      [--batch_col sample] [--gene_col ...] [--filter_col pass_guide_filter]
then:
  python scripts/17_xatlas_pseudobulk.py --h5ad <subset.h5ad> --out ... (same flags)
"""

import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

_spec = importlib.util.spec_from_file_location("pb17", Path(__file__).with_name("17_xatlas_pseudobulk.py"))
pb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pb)

ROWS_PER_CHUNK = 8000
BLOCK = 64 * 2 ** 20
MAX_RETRIES = 10


def open_remote(path):
    import h5py
    if "://" in path:
        import fsspec
        f = fsspec.open(path, "rb", block_size=BLOCK, cache_type="readahead").open()
        return h5py.File(f, "r")
    return h5py.File(path, "r")


def select_rows(h, a):
    obs, var = h["obs"], h["var"]
    target = pb.read_column(obs, a.target_col).astype(str)
    keep = np.ones(len(target), dtype=bool)
    if a.filter_col:
        keep &= pb.read_column(obs, a.filter_col).astype(bool)
    rng = np.random.default_rng(pb.SEED)                      # identical to script 17
    ctrl_all = np.where((target == a.control_label) & keep)[0]
    ctrl = np.sort(rng.choice(ctrl_all, min(pb.N_CTRL, len(ctrl_all)), replace=False))
    labels_avail = set(np.unique(target[keep]))
    labs = {pb.resolve(s, labels_avail) for s in pb.analysis_knockdowns()} - {None}
    kd_rows = np.where(np.isin(target, sorted(labs)) & keep)[0]
    rows = np.union1d(ctrl, kd_rows)
    batch = pb.read_column(obs, a.batch_col).astype(str) if a.batch_col else None
    genes = (pb.read_column(var, a.gene_col) if a.gene_col else pb.index_column(var)).astype(str)
    return rows, target, batch, genes, len(ctrl_all), len(ctrl)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--target_col", default="gene_target")
    p.add_argument("--control_label", required=True)
    p.add_argument("--batch_col")
    p.add_argument("--gene_col")
    p.add_argument("--filter_col")
    a = p.parse_args()
    t0 = time.time()
    work = Path(a.work)
    (work / "parts").mkdir(parents=True, exist_ok=True)

    h = open_remote(a.h5ad)
    X = h["X"]
    if "indptr" not in X:
        raise SystemExit("X is not CSR")
    n_rows, n_genes = (int(x) for x in X.attrs["shape"])
    sel_path = work / "selection.npz"
    if sel_path.exists():
        d = np.load(sel_path, allow_pickle=False)
        rows, target, genes = d["rows"], d["target"], d["genes"]
        batch = d["batch"] if "batch" in d.files else None
        n_ctrl_all, n_ctrl = int(d["n_ctrl_all"]), int(d["n_ctrl"])
    else:
        rows, target, batch, genes, n_ctrl_all, n_ctrl = select_rows(h, a)
        extra = {"batch": batch} if batch is not None else {}
        np.savez(sel_path, rows=rows, target=target, genes=genes, n_ctrl_all=n_ctrl_all, n_ctrl=n_ctrl, **extra)
    indptr = np.asarray(X["indptr"][()], dtype=np.int64)
    print(f"[extract] rows {n_rows}, genes {n_genes}, keep {len(rows)} "
          f"(controls {n_ctrl} of {n_ctrl_all}); {indptr[-1] / 1e9:.2f}e9 nonzeros", flush=True)

    starts = list(range(0, n_rows, ROWS_PER_CHUNK))
    for k, r0 in enumerate(starts):
        part = work / "parts" / f"part_{k:06d}.npz"
        if part.exists():
            continue
        r1 = min(r0 + ROWS_PER_CHUNK, n_rows)
        want = rows[(rows >= r0) & (rows < r1)]
        if len(want) == 0:
            np.savez(part, rows=want, lens=np.zeros(0, np.int64), data=np.zeros(0, np.float32),
                     indices=np.zeros(0, np.int32))
            continue
        for attempt in range(MAX_RETRIES):
            try:
                p0, p1 = indptr[r0], indptr[r1]
                data = np.asarray(X["data"][p0:p1], dtype=np.float32)
                ind = np.asarray(X["indices"][p0:p1], dtype=np.int32)
                break
            except Exception as e:  # network: back off, reopen, retry the same chunk
                wait = min(2 ** attempt, 120)
                print(f"[extract] chunk {k} attempt {attempt + 1} failed ({type(e).__name__}); retry in {wait}s",
                      flush=True)
                time.sleep(wait)
                h = open_remote(a.h5ad)
                X = h["X"]
        else:
            raise SystemExit(f"chunk {k} failed {MAX_RETRIES} times")
        lens, dd, ii = [], [], []
        for r in want:
            s, e = indptr[r] - p0, indptr[r + 1] - p0
            lens.append(e - s)
            dd.append(data[s:e])
            ii.append(ind[s:e])
        np.savez(part, rows=want, lens=np.array(lens, np.int64), data=np.concatenate(dd), indices=np.concatenate(ii))
        if (k + 1) % 20 == 0:
            done = (k + 1) / len(starts)
            el = time.time() - t0
            print(f"[extract] chunk {k + 1}/{len(starts)} ({100 * done:.1f}%), {el / 60:.0f} min, "
                  f"ETA {el / done * (1 - done) / 60:.0f} min", flush=True)

    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    R, L, D, I = [], [], [], []
    for k in range(len(starts)):
        d = np.load(work / "parts" / f"part_{k:06d}.npz")
        R.append(d["rows"]); L.append(d["lens"]); D.append(d["data"]); I.append(d["indices"])
    R, L = np.concatenate(R), np.concatenate(L)
    assert np.array_equal(R, rows), "extracted rows differ from the selection"
    M = sp.csr_matrix((np.concatenate(D), np.concatenate(I), np.concatenate([[0], np.cumsum(L)])),
                      shape=(len(R), n_genes))
    obs = {a.target_col: pd.Categorical(target[R])}
    if batch is not None:
        obs[a.batch_col] = pd.Categorical(batch[R])
    if a.filter_col:
        obs[a.filter_col] = np.ones(len(R), dtype=bool)        # only kept cells were extracted
    sub = ad.AnnData(X=M, obs=pd.DataFrame(obs, index=[f"r{r}" for r in R]),
                     var=pd.DataFrame(index=pd.Index(genes).astype(str)))
    sub.var_names_make_unique()
    sub.write_h5ad(a.out)
    meta = {"source": a.h5ad, "n_rows_source": n_rows, "n_rows_subset": int(len(R)),
            "n_controls_source": n_ctrl_all, "n_controls_subset": n_ctrl, "minutes": round((time.time() - t0) / 60, 1)}
    Path(a.out + ".meta.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
