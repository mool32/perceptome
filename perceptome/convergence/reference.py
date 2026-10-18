"""Bundled reference of the convergence paper: definitions, published results, SI tables.

Everything here is read from the manuscript (Spiro 2026, Zenodo 10.5281/zenodo.20542130):
the cell types that define the anchor / core / neighbourhood, the quality-control rules,
the cohort-level within-patient results, the retractions, and Supporting-Information
Tables S1, S2, S4, S6, S7 and S8 (as DataFrames). Nothing is recomputed from cancer data,
which are not part of this package.
"""

import json
from functools import lru_cache
from pathlib import Path

import pandas as pd

_DATA = Path(__file__).parent / "data"
_FILE = _DATA / "convergence_v1.json"

# index column of each bundled table
_TABLE_INDEX = {"S1": "cancer", "S2": "cell_type", "S4": "module", "S6": "cell_type",
                "S7": "module", "S8": "cancer"}


@lru_cache(maxsize=1)
def load_convergence_reference():
    """Load the bundled convergence reference.

    Returns
    -------
    dict
        source        citation of the manuscript and the bundled table files
        definitions   anchor / secondary, region definitions, quarantine, quality-control rules
        headline      contraction, enrichment, robustness and catalog-dependence numbers
        within_patient   list of cohort-level results (see ``published_within_patient``)
        retracted     statements the manuscript withdraws
        limits        what the manuscript states the data do not support
        tables        {'S1', 'S2', 'S4', 'S6', 'S7', 'S8'} DataFrames, indexed by cancer / cell type / module
    """
    with open(_FILE, encoding="utf-8") as f:
        ref = json.load(f)
    ref["tables"] = {key: pd.read_csv(_DATA / fname).set_index(_TABLE_INDEX[key])
                     for key, fname in ref["source"]["tables"].items()}
    return ref


def published_within_patient():
    """Cohort-level within-patient results of the paper as a flat DataFrame.

    Compare your own ``within_patient_occupancy`` output with these. Shape-only and
    complexity-matched variants are in separate columns where the paper reports them.
    """
    rows = []
    for rec in load_convergence_reference()["within_patient"]:
        shape = rec.get("shape_only") or {}
        rows.append({
            "id": rec["id"], "tissue": rec["tissue"], "cohort": rec["cohort"], "design": rec["design"],
            "n_patients": rec.get("n_patients", rec.get("n_donors")),
            "n_deeper": rec.get("n_deeper"),
            "normal_pct": rec.get("normal_pct", rec.get("control_pct")),
            "malignant_pct": rec.get("malignant_pct", rec.get("disease_pct")),
            "statistic": rec.get("statistic"),
            "shape_only_delta_points": shape.get("delta_points"),
            "verdict": rec.get("verdict"),
        })
    return pd.DataFrame(rows).set_index("id")


def region_members(region="core", ranking=None):
    """Cell types of a convergence region.

    Parameters
    ----------
    region : 'anchor' | 'core' | 'three_cell_core' | 'neighbourhood' | 'region' | int | sequence of str
        ``anchor`` (cytotrophoblast), ``core`` (cytotrophoblast + megakaryocyte),
        ``three_cell_core`` (+ migrating cytotrophoblast), ``neighbourhood`` (the four-cell
        majority-enriched set) are fixed lists. ``region`` / an integer ``k`` take the top ``k``
        candidates of the cross-cancer recurrence ranking (default k = 14; ties at the boundary
        are broken by median then mean depth — Table S2 does not resolve places 13-14).
        A sequence of names is validated against the reference and returned as a list.
    ranking : DataFrame | None
        A recurrence ranking from :func:`recurrence_ranking`; default is the paper's Table S2.
    """
    ref = load_convergence_reference()
    defs = ref["definitions"]
    if isinstance(region, str) and region in defs["regions"]:
        return list(defs["regions"][region])
    if isinstance(region, str) and region == "region":
        region = defs["region_rule"]["k"]
    if isinstance(region, (int,)) and not isinstance(region, bool):
        rk = ref["tables"]["S2"] if ranking is None else ranking
        rk = rk[rk["qc_status"] == "candidate"] if "qc_status" in rk.columns else rk
        by = [c for c in defs["region_rule"]["rank_by"] if c in rk.columns]
        top = rk.sort_values(by, ascending=False, kind="mergesort").head(region)
        return list(top.index)
    if isinstance(region, str):
        raise ValueError(f"Unknown region {region!r}; use one of {sorted(defs['regions'])}, 'region', an int, or a list of cell types.")
    members = list(region)
    known = set(ref["tables"]["S2"].index)
    unknown = [m for m in members if m not in known]
    if unknown:
        raise ValueError(f"Not reference cell types: {unknown}")
    return members


def quarantined_cell_types():
    """The two post-cellular low-complexity reference cells removed before assignment."""
    return list(load_convergence_reference()["definitions"]["quarantine"])
