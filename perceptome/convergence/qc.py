"""Cell-level quality control of the paper's cancer projection (Methods, Note S1, Note S3).

Per cohort: uppercase and de-duplicate gene names; keep cells with at least 1,000 detected genes;
``normalize_total(1e4)`` + ``log1p`` if the matrix is not already log-normalized (max > 100).
The 1,000-gene floor exists because module scores of low-complexity cells are degenerate and
map preferentially to low-complexity reference types (Note S3).
"""

import numpy as np


def detected_genes(X):
    """Number of genes with non-zero expression per cell (sparse or dense)."""
    return np.asarray((X > 0).sum(axis=1)).ravel()


def prepare_expression(adata, min_genes=1000, normalize=True, uppercase_genes=True,
                       target_sum=1e4, log_threshold=100, copy=True):
    """Apply the paper's per-cell preprocessing.

    Parameters
    ----------
    adata : AnnData
        Cells x genes, raw counts or log-normalized. Genes must be HGNC symbols (map Ensembl
        identifiers first).
    min_genes : int
        Minimum detected genes per cell (counted on the matrix as given, before normalization).
    normalize : bool
        If True and ``X.max() > log_threshold`` the matrix is treated as counts and
        ``normalize_total(target_sum)`` + ``log1p`` are applied.
    uppercase_genes : bool
        Uppercase gene symbols and keep the first of any duplicates. Numeric suffixes are *not*
        stripped: ``NKX3-1`` is a real symbol (a loader that undid ``var_names_make_unique`` by
        stripping ``-<digits>`` silently dropped it from every cohort, Note S-Q).

    Returns
    -------
    (AnnData, dict)
        The filtered, normalized object and a report with ``n_input``, ``n_pass``, ``min_genes``,
        ``median_genes_input``, ``normalized``.
    """
    import scanpy as sc

    ad = adata.copy() if copy else adata
    if uppercase_genes:
        ad.var_names = [str(g).upper() for g in ad.var_names]
        ad = ad[:, ~ad.var_names.duplicated(keep="first")].copy()

    n_genes = detected_genes(ad.X)
    keep = n_genes >= min_genes
    report = {"n_input": int(ad.n_obs), "n_pass": int(keep.sum()), "min_genes": int(min_genes),
              "median_genes_input": float(np.median(n_genes)) if len(n_genes) else float("nan")}
    ad = ad[keep].copy()
    ad.obs["n_genes_detected"] = n_genes[keep]

    xmax = ad.X.max() if ad.n_obs else 0
    report["normalized"] = bool(normalize and xmax > log_threshold)
    if report["normalized"]:
        sc.pp.normalize_total(ad, target_sum=target_sum)
        sc.pp.log1p(ad)
    return ad, report
