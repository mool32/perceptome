"""Assign cells to their nearest normal cell type and measure occupancy of a region.

The paper's convergence statistics all rest on one operation (Methods, "Cancer projection"):
a cell is scored on the 44 modules, standardized against the reference calibration, projected
through the fixed loadings, and assigned to the nearest of the normal cell-type centroids by
cosine. The erythrocyte and the platelet — post-cellular, low-complexity bodies — are removed
from the candidate pool *before* that step; the cells that would have gone to them are
reassigned to their best surviving neighbour (Note S-K), leaving 152 candidates.
"""

import numpy as np
import pandas as pd

from ..eigenspace.calibration import Calibration, load_calibration, warn_if_scale_mismatch
from .reference import load_convergence_reference, region_members, quarantined_cell_types


def _unit_rows(X):
    norm = np.linalg.norm(X, axis=1, keepdims=True)
    safe = np.where(norm < 1e-12, np.nan, norm)
    return X / safe


def assign_to_reference(scores, calibration="pseudobulk", shape_only=False, quarantine=True,
                        missing_modules=None, warn_scale=True, return_coordinates=False):
    """Top-1 nearest normal cell type (cosine) for each row of module scores.

    Parameters
    ----------
    scores : DataFrame
        cells x modules (``mean_raw``, from ``pct.score_modules``).
    calibration : 'pseudobulk' | 'single_cell_scaled' | Calibration
        Reference scaling. **Per-cell data need ``single_cell_scaled``**: under the pseudobulk
        calibration cells collapse onto low-complexity targets (Notes S-B, S-H). The pseudobulk
        calibration is bundled; the single-cell one has to be built from the HPA table
        (``scripts/13_build_single_cell_calibration.py``).
    shape_only : bool
        Remove each cell's overall module magnitude first (subtract its own mean over the 44
        modules, reference and cells alike; Note S-K iv).
    quarantine : bool
        Remove the erythrocyte and platelet from the candidate pool (default, as in the paper).
    missing_modules : iterable of str | None
        Modules to treat as unmeasured (placed at the reference mean).
    return_coordinates : bool
        Also return the standardized eigenspace coordinates of the cells.

    Returns
    -------
    DataFrame  (or (DataFrame, coordinates) with ``return_coordinates``)
        index = cells; columns ``nearest`` (reference cell type), ``cosine`` (to it),
        ``first_choice`` (nearest among all 154, before quarantine), ``displaced`` (True when
        the first choice was a quarantined cell type and the cell was reassigned).
    """
    cal = calibration if isinstance(calibration, Calibration) else load_calibration(calibration)
    if shape_only:
        cal = cal.shape_only_calibration()
    coords = cal.project(scores, missing_modules=missing_modules)
    if warn_scale:
        warn_if_scale_mismatch(cal.standardize(scores), cal, stacklevel=3)

    centroids = cal.centroids
    sim_all = _unit_rows(coords.values) @ _unit_rows(centroids.values).T
    names = np.array(centroids.index)
    quar = set(quarantined_cell_types()) if quarantine else set()
    cand = np.array([n not in quar for n in names])

    sim = np.where(cand[None, :], sim_all, -np.inf)
    ok = ~np.isnan(sim_all).all(axis=1)
    best = np.where(ok, np.nan_to_num(sim, nan=-np.inf).argmax(axis=1), 0)
    first = np.where(ok, np.nan_to_num(sim_all, nan=-np.inf).argmax(axis=1), 0)

    out = pd.DataFrame({
        "nearest": np.where(ok, names[best], None),
        "cosine": np.where(ok, sim_all[np.arange(len(coords)), best], np.nan),
        "first_choice": np.where(ok, names[first], None),
    }, index=scores.index)
    out["displaced"] = out["first_choice"].isin(quar)
    out.attrs["calibration"] = cal.name
    return (out, coords) if return_coordinates else out


def assign_cells(adata, calibration="pseudobulk", min_genes=1000, shape_only=False,
                 normalize=True, resolve_aliases=False, prepare=True):
    """The paper's per-cell pipeline: QC, score 44 modules, standardize, project, assign.

    Parameters
    ----------
    adata : AnnData
        Cells x genes (counts or log-normalized), HGNC symbols.
    prepare : bool
        Apply :func:`perceptome.convergence.prepare_expression` (uppercase genes, >= ``min_genes``
        detected genes, normalize_total(1e4)+log1p if the matrix holds counts).
    resolve_aliases : bool
        Forwarded to ``score_modules``: also match current symbols of the renamed catalog genes
        (MB21D1/TMEM173/ARNTL -> CGAS/STING1/BMAL1). Off by default because the bundled reference
        was built without it (see docs/PAPER_MAP.md).

    Returns
    -------
    dict
        assignment    DataFrame from :func:`assign_to_reference` (cells that passed QC)
        scores        module scores (cells x 44)
        coordinates   standardized eigenspace coordinates
        coverage      per-module gene coverage
        qc            QC report
    """
    from ..score import score_modules

    if prepare:
        from .qc import prepare_expression
        adata, qc = prepare_expression(adata, min_genes=min_genes, normalize=normalize)
    else:
        qc = {"n_input": int(adata.n_obs), "n_pass": int(adata.n_obs), "min_genes": None, "normalized": False}

    scored = score_modules(adata, method="mean_raw", resolve_aliases=resolve_aliases)
    scores, coverage = scored["scores"], scored["coverage"]
    uncovered = list(coverage.index[coverage["n_genes_found"] < 2])

    cal = calibration if isinstance(calibration, Calibration) else load_calibration(calibration)
    assignment, coords = assign_to_reference(scores, calibration=cal, shape_only=shape_only,
                                             missing_modules=uncovered, return_coordinates=True)
    return {"assignment": assignment, "scores": scores, "coordinates": coords,
            "coverage": coverage, "qc": qc}


def _nearest(assignment):
    return assignment["nearest"] if isinstance(assignment, pd.DataFrame) else pd.Series(assignment)


def in_region(assignment, region="core", ranking=None):
    """Boolean Series: is each cell's nearest cell type a member of ``region``?"""
    members = set(region_members(region, ranking=ranking))
    return _nearest(assignment).isin(members)


def core_fraction(assignment, region="core", ranking=None):
    """Fraction of cells whose nearest normal cell type lies in ``region`` (default: the two-cell core)."""
    hit = in_region(assignment, region, ranking=ranking)
    return float(hit.mean()) if len(hit) else float("nan")


def enrichment(fraction, k, n_candidates=None):
    """Fold enrichment of an occupancy fraction over the ``k / n_candidates`` chance baseline."""
    if n_candidates is None:
        n_candidates = load_convergence_reference()["definitions"]["n_candidates"]
    return float(fraction) / (k / n_candidates)


def region_occupancy(assignment, regions=("anchor", "core", "neighbourhood", "region")):
    """Occupancy and enrichment of several boundaries (the paper's Fig. 3F for your cells)."""
    rows = []
    for name in regions:
        members = region_members(name)
        frac = core_fraction(assignment, name)
        rows.append({"region": name, "k": len(members), "occupancy_pct": 100 * frac,
                     "chance_pct": 100 * len(members) / load_convergence_reference()["definitions"]["n_candidates"],
                     "fold_over_chance": enrichment(frac, len(members))})
    return pd.DataFrame(rows).set_index("region")
