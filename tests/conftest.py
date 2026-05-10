"""Shared test fixtures."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


# Make in-tree perceptome importable
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


@pytest.fixture
def tiny_adata():
    """Synthetic AnnData with 50 cells × 60 genes covering several modules.

    Includes core/activity genes for HSF1, UPR-ATF6, NPAS4, ERK/MAPK, Cell Cycle.
    cells split into two cell types (A: 30, B: 20) and two conditions (ctrl: 25, drug: 25).
    """
    try:
        import anndata
    except ImportError:
        pytest.skip("anndata not installed")

    rng = np.random.default_rng(0)

    genes = [
        # HSF1 core + activity
        "HSF1", "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "DNAJB1", "BAG3",
        # UPR-ATF6 core + activity
        "ATF6", "HSPA5", "HSP90B1", "CALR", "PDIA4", "XBP1",
        # NPAS4 core + activity
        "NPAS4", "ARNT2", "CAMK2A", "CAMK2B", "BDNF", "NPTX2", "GABRA1",
        # ERK/MAPK core + activity
        "MAPK1", "MAPK3", "FOS", "EGR1", "DUSP1", "DUSP6", "MYC",
        # Cell Cycle (Tirosh)
        "MKI67", "PCNA", "TOP2A", "MCM2", "CCNA2", "CCNB1", "CDK1", "AURKA",
        # Housekeeping
        "ACTB", "GAPDH", "PPIA", "RPL13A", "RPLP0", "TBP", "HPRT1", "PGK1",
        # Filler
        *[f"GENE{i}" for i in range(60 - 7 - 6 - 7 - 7 - 8 - 8)],
    ]
    n_cells = 50
    n_genes = len(genes)
    X = rng.normal(loc=2.0, scale=0.5, size=(n_cells, n_genes))
    X = np.clip(X, 0, None)

    obs = pd.DataFrame({
        "cell_type": ["A"] * 30 + ["B"] * 20,
        "tissue": ["epithelial"] * 30 + ["immune"] * 20,
        "condition": ["ctrl"] * 25 + ["drug"] * 25,
        "drug_name": ["DMSO"] * 25 + ["trametinib"] * 25,
    }, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=genes)
    return anndata.AnnData(X=X, obs=obs, var=var)
