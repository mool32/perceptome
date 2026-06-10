"""Varimax (orthogonal) rotation of the eigenspace loadings.

Coordinate-invariance robustness, lens #1 (lowest risk). Varimax is an ORTHOGONAL
rotation: it preserves every pairwise distance in the full PC space, so all
distance-based findings (cluster compactness, nearest-reference, cosine of shift
vectors) are *automatically invariant* to it. It only redistributes variance so
each module loads on fewer axes — cleaner axes for interpretation and 2-D views,
not a different geometry. PCA stays the canonical space.
"""

import numpy as np


def varimax(loadings, gamma=1.0, max_iter=200, tol=1e-7):
    """Kaiser varimax rotation.

    Parameters
    ----------
    loadings : ndarray (n_modules × n_pcs)
        Eigenspace loadings (e.g. from reference_v03 loadings matrix).
    gamma : float
        1.0 = varimax. (gamma=0 → quartimax.)
    max_iter, tol : int, float
        Convergence controls.

    Returns
    -------
    rotated : ndarray (n_modules × n_pcs)   rotated loadings
    R       : ndarray (n_pcs × n_pcs)        orthogonal rotation matrix (rotated = loadings @ R)
    """
    L = np.asarray(loadings, dtype=float)
    p, k = L.shape
    if k < 2:
        return L.copy(), np.eye(k)
    R = np.eye(k)
    d = 0.0
    for _ in range(max_iter):
        Lam = L @ R
        diag = np.diag((Lam ** 2).sum(axis=0))
        u, s, vt = np.linalg.svd(L.T @ (Lam ** 3 - (gamma / p) * Lam @ diag))
        R = u @ vt
        d_old = d
        d = s.sum()
        if d_old != 0 and abs(d - d_old) / d_old < tol:
            break
    return L @ R, R
