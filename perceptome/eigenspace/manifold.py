"""Diffusion-map embedding of the perceptome reference manifold.

Coordinate-invariance robustness, lens #2 (the real nonlinear test). PCA sees only
linear geometry; if the 154 normal cell states lie on a curved manifold (as
developmental trajectories usually do), a region like the cancer-attractor beacon
cluster may be sharper on the manifold than in any PCA plane.

Why diffusion maps and not UMAP/t-SNE (per the project's geometry note):
  * diffusion distance is a TRUE metric — nearest-reference and compactness
    queries are meaningful (UMAP/t-SNE preserve only local topology);
  * out-of-sample projection is principled (Nyström extension) — cancer cells
    enter by projection into a space fit on normal cells only, never refit;
  * bandwidth ε is fixed by a DECLARED rule (median pairwise sq-distance), not
    tuned — so there is no "turn the knob until the cluster appears" risk.

NB: diffusion maps were already shown to lose to PCA as a *classifier* (theory
overview §III: kNN 0.494 vs 0.552). This module is NOT a replacement — it is a
reproducibility lens: does a finding seen in PCA survive in an independent,
deliberately weaker geometry?
"""

import numpy as np


def _sqdist(A, B=None):
    B = A if B is None else B
    a2 = (A ** 2).sum(1)[:, None]
    b2 = (B ** 2).sum(1)[None, :]
    return np.maximum(a2 + b2 - 2.0 * A @ B.T, 0.0)


class DiffusionMap:
    """Coifman–Lafon diffusion map. Fit on a reference set; project new points via Nyström.

    Parameters
    ----------
    n_components : int
        Non-trivial diffusion coordinates to keep.
    alpha : float
        Density normalization exponent (1.0 = recover the Laplace–Beltrami
        operator, removing the effect of uneven cell-type sampling; 0.0 = graph
        Laplacian).
    t : int
        Diffusion time (coords scaled by λ**t).
    epsilon : "median" | float
        Kernel bandwidth. "median" = median pairwise squared distance of the fit
        set (declared, not tuned).
    """

    def __init__(self, n_components=8, alpha=1.0, t=1, epsilon="median"):
        self.n_components = n_components
        self.alpha = alpha
        self.t = t
        self.epsilon = epsilon

    def fit(self, X):
        X = np.asarray(X, dtype=float)
        n = len(X)
        D2 = _sqdist(X)
        if self.epsilon == "median":
            eps = float(np.median(D2[np.triu_indices(n, 1)]))
        else:
            eps = float(self.epsilon)
        eps = eps if eps > 0 else 1.0
        K = np.exp(-D2 / eps)
        q = K.sum(1)
        Ka = K / np.outer(q ** self.alpha, q ** self.alpha)   # α-normalization
        d = Ka.sum(1)
        dinv_sqrt = 1.0 / np.sqrt(d)
        Ms = Ka * np.outer(dinv_sqrt, dinv_sqrt)               # symmetric conjugate of P
        evals, evecs = np.linalg.eigh(Ms)
        order = np.argsort(evals)[::-1]
        evals, evecs = evals[order], evecs[:, order]
        psi = evecs * dinv_sqrt[:, None]                       # right eigenvectors of P

        self.X_fit_, self.epsilon_ = X, eps
        self.q_fit_, self.evals_, self.psi_ = q, evals, psi
        kc = self.n_components
        self.coords_ = psi[:, 1:kc + 1] * (evals[1:kc + 1] ** self.t)
        return self

    def transform(self):
        """Diffusion coordinates of the fit set (n_fit × n_components)."""
        return self.coords_

    def transform_oos(self, X_new):
        """Nyström out-of-sample projection of new points into the fit manifold."""
        X_new = np.asarray(X_new, dtype=float)
        K = np.exp(-_sqdist(X_new, self.X_fit_) / self.epsilon_)
        qx = K.sum(1)
        Ka = K / np.outer(qx ** self.alpha, self.q_fit_ ** self.alpha)
        P = Ka / Ka.sum(1)[:, None]                            # n_new × n_fit transition rows
        kc = self.n_components
        psi_new = np.column_stack([
            (P @ self.psi_[:, l]) / self.evals_[l] for l in range(1, kc + 1)
        ])
        return psi_new * (self.evals_[1:kc + 1] ** self.t)
