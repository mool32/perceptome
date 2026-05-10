"""Build a new reference eigenspace from a cell-type × module activity matrix.

Used by scripts/03_build_eigenspace.py to produce the bundled v0.3 reference
from the 154 × 44 HPA matrix. Also exposed as `pct.rebuild()` for users who
want to derive a tissue- or condition-specific eigenspace from their own data.
"""

import numpy as np


def rebuild(module_scores_matrix, kaiser=True, bootstrap_n=100, seed=42):
    """Compute eigendecomposition + bootstrap stability for a cell-type × module matrix.

    Parameters
    ----------
    module_scores_matrix : DataFrame
        Cell types × modules (e.g., 154 × 44).
    kaiser : bool
        If True, n_kaiser = number of eigenvalues > 1.
    bootstrap_n : int
        Number of bootstrap resamples for per-PC stability.
    seed : int
        RNG seed for reproducibility.

    Returns
    -------
    dict   serializable, ready to be saved as reference_v03.json
        eigenvalues          list[float]
        n_kaiser             int
        stability            dict {PC: bootstrap-cosine}
        spectral_fit         {alpha, R2}    log-log power-law fit
        module_order         list[str]
        loadings             {PCi: {module: loading}}
        var_explained        list[float]
        n_pcs                int
        n_cell_types_used    int
        version              str
    """
    rng = np.random.default_rng(seed)

    X = module_scores_matrix.values.astype(float)
    n_ct, n_mod = X.shape

    means = X.mean(axis=0)
    stds = X.std(axis=0)
    stds[stds < 1e-10] = 1.0
    Z = (X - means) / stds

    corr = np.corrcoef(Z.T)
    eigvals, eigvecs = np.linalg.eigh(corr)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    n_kaiser = int(np.sum(eigvals > 1.0)) if kaiser else n_mod

    stability = {}
    if bootstrap_n > 0:
        boot = np.zeros((bootstrap_n, n_mod, n_mod))
        for b in range(bootstrap_n):
            ix = rng.choice(n_ct, size=n_ct, replace=True)
            Zb = Z[ix]
            cb = np.corrcoef(Zb.T)
            ev, vec = np.linalg.eigh(cb)
            sort = np.argsort(ev)[::-1]
            boot[b] = vec[:, sort]
        for pc in range(min(n_kaiser, n_mod)):
            ref_vec = eigvecs[:, pc]
            cosines = []
            for b in range(bootstrap_n):
                bv = boot[b, :, pc]
                cos = abs(np.dot(ref_vec, bv) / (np.linalg.norm(ref_vec) * np.linalg.norm(bv) + 1e-12))
                cosines.append(cos)
            stability[f"PC{pc+1}"] = float(np.mean(cosines))

    from scipy import stats as sp_stats
    k = np.arange(1, len(eigvals) + 1)
    pos = eigvals > 0
    if pos.sum() > 2:
        slope, intercept, r, _, _ = sp_stats.linregress(np.log(k[pos]), np.log(eigvals[pos]))
        spectral_fit = {"alpha": float(-slope), "R2": float(r ** 2)}
    else:
        spectral_fit = {"alpha": None, "R2": None}

    module_order = list(module_scores_matrix.columns)
    loadings = {
        f"PC{pc+1}": {module_order[i]: float(eigvecs[i, pc]) for i in range(n_mod)}
        for pc in range(n_kaiser)
    }
    var_explained = (eigvals / eigvals.sum()).tolist()

    return {
        "eigenvalues": eigvals.tolist(),
        "n_kaiser": n_kaiser,
        "stability": stability,
        "bootstrap_stability": [stability.get(f"PC{i+1}", 0.0) for i in range(n_kaiser)],
        "spectral_fit": spectral_fit,
        "module_order": module_order,
        "loadings": loadings,
        "var_explained": var_explained,
        "n_pcs": n_kaiser,
        "n_cell_types_used": n_ct,
        "n_modules": n_mod,
        "mean_per_module": means.tolist(),
        "std_per_module": stds.tolist(),
    }
