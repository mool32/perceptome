"""Compare sample positions to reference vector databases.

References supported (under perceptome/reference/data/):
  - disease   AD/DKD/RA/IPF perturbation vectors
  - aging     inflammaging + collapse axes
  - drugs     CMap/LINCS L1000 (recomputed for v0.3 in last commit)
  - attractor Cancer capacity-direction (Paper 4.2 P3, NEW in v0.3)
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..utils import cosine

_REF_DIR = Path(__file__).resolve().parent.parent / "reference" / "data"


def compare_to_references(
    coordinates,
    include_diseases=True,
    include_aging=True,
    include_attractor=True,
    references_dir=None,
):
    """Compute cosine similarities to reference perturbation vectors.

    NOTE: Drug-disease cosine matching was removed in v0.2 — Paper 4.1 falsified
    it in three independent formulations (cosine v1.0, subspace v1.1, module
    overlap v1.2). The validated drug operation is activity-layer screen against
    a background null on TF target panels — use pct.activity_layer_screen and
    pct.drug_anchors instead.

    Parameters
    ----------
    coordinates : DataFrame
        Samples × PCs (output of project()).
    include_diseases / include_aging / include_attractor : bool
    references_dir : str | Path | None
        Override the bundled reference directory.

    Returns
    -------
    dict
        disease_similarities    DataFrame (samples × diseases)
        cos_inflammaging        Series
        cos_collapse            Series
        aging_type              Series  (inflammaging|perception_collapse|mixed|neutral)
        cos_attractor           Series  (cosine with cancer capacity-direction)
        attractor_quadrant      Series  (high|medium|low alignment)
        eigenspace_region       Series  categorical PC1/PC2 region label
    """
    ref_dir = Path(references_dir) if references_dir else _REF_DIR
    n_pcs = coordinates.shape[1]
    out = {}

    if include_diseases:
        out.update(_disease_similarities(coordinates, ref_dir, n_pcs))

    if include_aging:
        out.update(_aging_similarities(coordinates, ref_dir, n_pcs))

    if include_attractor:
        out.update(_attractor_alignment(coordinates, ref_dir, n_pcs))

    out["eigenspace_region"] = _classify_region(coordinates)
    return out


def _disease_similarities(coordinates, ref_dir, n_pcs):
    p = ref_dir / "disease_vectors.json"
    if not p.exists():
        # v0.2 vectors are pinned as disease_vectors_v02pinned.json (12-PC v0.2 space,
        # incompatible with 9-PC v0.3 eigenspace). Recompute scheduled for v0.2.1.
        return {"disease_similarities_status": "pending_v0.2.1_recompute"}
    with open(p) as f:
        diseases = json.load(f)
    sims = {}
    for name, d in diseases.items():
        v = d.get("mean_vector_12D") or d.get("mean_vector")
        if v is None:
            continue
        v = np.array(v[:n_pcs])
        sims[name] = [cosine(row.values[:n_pcs], v) for _, row in coordinates.iterrows()]
    return {"disease_similarities": pd.DataFrame(sims, index=coordinates.index)}


def _aging_similarities(coordinates, ref_dir, n_pcs):
    p = ref_dir / "aging_reference.json"
    if not p.exists():
        # v0.2 aging vectors pinned as aging_reference_v02pinned.json (12-PC v0.2 space).
        # Recompute scheduled for v0.2.1.
        return {"aging_similarities_status": "pending_v0.2.1_recompute"}
    with open(p) as f:
        aging_ref = json.load(f)
    res = {}
    inf_v = aging_ref.get("inflammaging_direction")
    col_v = aging_ref.get("collapse_direction")
    if inf_v:
        v = np.array(inf_v[:n_pcs])
        res["cos_inflammaging"] = pd.Series(
            [cosine(r.values[:n_pcs], v) for _, r in coordinates.iterrows()],
            index=coordinates.index, name="cos_inflammaging",
        )
    if col_v:
        v = np.array(col_v[:n_pcs])
        res["cos_collapse"] = pd.Series(
            [cosine(r.values[:n_pcs], v) for _, r in coordinates.iterrows()],
            index=coordinates.index, name="cos_collapse",
        )
    if "cos_inflammaging" in res and "cos_collapse" in res:
        labels = []
        for i in range(len(coordinates)):
            inf = res["cos_inflammaging"].iloc[i]
            col = res["cos_collapse"].iloc[i]
            if inf > 0.3 and col < 0:
                labels.append("inflammaging")
            elif col > 0.3 and inf < 0:
                labels.append("perception_collapse")
            elif inf > 0.2 and col > 0.2:
                labels.append("mixed")
            else:
                labels.append("neutral")
        res["aging_type"] = pd.Series(labels, index=coordinates.index, name="aging_type")
    return res


def _attractor_alignment(coordinates, ref_dir, n_pcs):
    p = ref_dir / "attractor_v1.json"
    if not p.exists():
        return {}
    with open(p) as f:
        attr = json.load(f)
    v = attr.get("attractor_direction_eigenspace")
    if not v:
        return {}
    v = np.array(v[:n_pcs])
    cos = pd.Series(
        [cosine(r.values[:n_pcs], v) for _, r in coordinates.iterrows()],
        index=coordinates.index, name="cos_attractor",
    )
    quad = cos.apply(
        lambda x: "high" if x > 0.3 else "low" if x < -0.3 else "medium"
    ).rename("attractor_quadrant")
    return {"cos_attractor": cos, "attractor_quadrant": quad}


def _classify_region(coordinates):
    """Coarse PC1/PC2 region label.

    NOTE: Cancer detection requires tissue-matched comparison via compare_conditions(),
    not absolute PC1 — cross-tissue PC1 differences (~3-5) dwarf cancer signal (~0.06).
    """
    labels = []
    for _, r in coordinates.iterrows():
        pc1 = r.get("PC1", 0)
        pc2 = r.get("PC2", 0)
        if pc1 > 2:
            labels.append("high_perception")
        elif pc1 < -1:
            labels.append("low_perception")
        elif pc2 > 1.5:
            labels.append("stress_response")
        elif pc2 < -1.5:
            labels.append("growth_mode")
        else:
            labels.append("moderate")
    return pd.Series(labels, index=coordinates.index, name="region")
