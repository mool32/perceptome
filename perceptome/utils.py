"""Shared utilities (cosine, gene-set means, paths)."""

from pathlib import Path
import numpy as np


def cosine(u, v):
    """Cosine similarity between two 1-D arrays. Returns 0 for zero-norm."""
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu < 1e-12 or nv < 1e-12:
        return 0.0
    return float(np.dot(u, v) / (nu * nv))


def package_root() -> Path:
    """Return the perceptome package root directory."""
    return Path(__file__).parent
