"""Per-cell module scoring methods.

Three methods, by recommended use:
  mean_raw          (default)         absolute log-norm levels; safe across datasets
  scanpy_corrected  (proliferation)   background-corrected via scanpy.tl.score_genes;
                                      USE THIS when comparing cells with very different
                                      transcriptome complexity (e.g., proliferating vs
                                      post-mitotic). Caught by Paper 4.5 v1.2 amendment
                                      after random-200 ARTIFACT FAIL detection.
  mean_zscore       (within-dataset)  z across cells; DO NOT compare two independently
                                      z-scored datasets — destroys cross-dataset signal.

Two convenience wrappers expose the core/activity dichotomy explicitly:
  score_readiness  → score_modules(gene_set='core')   (R)
  score_activity   → score_modules(gene_set='activity') (A)

These are the inputs to pct.perceptivity.compute().
"""

import numpy as np
import pandas as pd

from ..catalog import load_catalog, get_genes


def score_modules(
    adata,
    catalog=None,
    method="mean_raw",
    gene_set="core",
    min_genes=2,
    return_per_cell=True,
    cluster_column=None,
):
    """Compute per-cell (or per-cluster) module activity scores.

    Parameters
    ----------
    adata : AnnData
        Cells × genes (log-normalized expression in .X expected).
    catalog : dict | None
        Module catalog. Default = bundled 44-module v0.3.
    method : 'mean_raw' | 'scanpy_corrected' | 'mean_zscore'
    gene_set : 'core' | 'activity'
        Which gene set to score (machinery presence vs TF target engagement).
    min_genes : int
        Minimum genes that must be found in adata.var_names per module.
    return_per_cell : bool
        If False and cluster_column given, return per-cluster means.
    cluster_column : str | None
        adata.obs column for cluster aggregation when return_per_cell=False.

    Returns
    -------
    dict
        scores         DataFrame (cells/clusters × n_modules)
        coverage       DataFrame (modules × {n_genes_total, n_genes_found, coverage, genes_found})
        missing_genes  dict {module: [missing genes]}
        cognitive_load Series (per cell/cluster: count of modules above the across-cells median)
    """
    cat = load_catalog() if (catalog is None or isinstance(catalog, str)) else catalog
    if isinstance(catalog, str):
        cat = load_catalog(catalog)

    modules = sorted(cat["modules"].keys())
    var_set = set(adata.var_names)

    scores_dict = {}
    cov_rows = []
    missing = {}

    for mod in modules:
        genes = get_genes(mod, gene_set=gene_set, catalog=cat)
        found = [g for g in genes if g in var_set]
        notf = [g for g in genes if g not in var_set]
        missing[mod] = notf
        cov_rows.append({
            "module": mod,
            "n_genes_total": len(genes),
            "n_genes_found": len(found),
            "coverage": len(found) / len(genes) if genes else 0.0,
            "genes_found": ",".join(found),
        })

        if len(found) < min_genes:
            scores_dict[mod] = np.zeros(adata.n_obs, dtype=float)
            continue

        if method == "scanpy_corrected":
            import scanpy as sc
            col = f"_pct_{mod}"
            try:
                sc.tl.score_genes(adata, gene_list=found, score_name=col, use_raw=False)
                scores_dict[mod] = adata.obs[col].values.astype(float)
                del adata.obs[col]
            except (RuntimeError, ValueError):
                # Fallback for tiny gene pools where score_genes refuses
                X_sub = _to_dense(adata[:, found].X)
                scores_dict[mod] = np.mean(X_sub, axis=1)
            continue

        X_sub = _to_dense(adata[:, found].X)
        cell_means = np.mean(X_sub, axis=1)

        if method == "mean_raw":
            scores_dict[mod] = cell_means
        elif method == "mean_zscore":
            mu = float(np.mean(cell_means))
            sd = float(np.std(cell_means))
            scores_dict[mod] = (cell_means - mu) / sd if sd > 1e-10 else np.zeros_like(cell_means)
        else:
            raise ValueError(
                f"method must be one of mean_raw|scanpy_corrected|mean_zscore, got {method!r}"
            )

    scores_df = pd.DataFrame(scores_dict, index=adata.obs_names)

    if not return_per_cell and cluster_column:
        scores_df = scores_df.groupby(adata.obs[cluster_column]).mean()

    medians = scores_df.median(axis=0)
    cog_load = (scores_df > medians).sum(axis=1)
    cog_load.name = "cognitive_load"

    return {
        "scores": scores_df,
        "coverage": pd.DataFrame(cov_rows).set_index("module"),
        "missing_genes": missing,
        "cognitive_load": cog_load,
    }


def score_readiness(adata, **kwargs):
    """Convenience: score_modules(gene_set='core'). Returns scores DataFrame only."""
    kwargs.setdefault("method", "mean_raw")
    return score_modules(adata, gene_set="core", **kwargs)["scores"]


def score_activity(adata, **kwargs):
    """Convenience: score_modules(gene_set='activity'). Returns scores DataFrame only."""
    kwargs.setdefault("method", "mean_raw")
    return score_modules(adata, gene_set="activity", **kwargs)["scores"]


def _to_dense(X):
    """Coerce sparse/dense matrix to a dense float ndarray."""
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=float)
