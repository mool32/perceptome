"""Compare two conditions (disease vs control, treated vs untreated)."""

import json
import numpy as np
import pandas as pd
from scipy import stats

from ..utils import cosine
from ..eigenspace.project import _load_reference
from .regime import infrastructure_regime


def compare_conditions(
    scores_condition1,
    scores_condition2,
    cell_type_column=None,
    obs1=None,
    obs2=None,
    project_to_eigenspace=True,
    reference=None,
):
    """Per-module + eigenspace delta between two conditions.

    Returns a dict with delta_modules (per-module test stats), delta_eigenspace
    (12-PC perturbation vector), magnitude, top_modules, disease_type, and
    infrastructure_regime. If cell_type_column / obs1 / obs2 provided, also
    returns per_cell_type breakdown plus cross-cell-type cosine convergence.
    """
    modules = sorted(scores_condition1.columns)

    def _compare_one(s1, s2):
        deltas = []
        for mod in modules:
            v1 = s1[mod].values.astype(float)
            v2 = s2[mod].values.astype(float)
            d = float(np.mean(v1) - np.mean(v2))
            try:
                _, pval = stats.mannwhitneyu(v1, v2, alternative="two-sided")
            except Exception:
                pval = 1.0
            pooled = float(np.sqrt((np.var(v1) + np.var(v2)) / 2))
            es = d / pooled if pooled > 1e-10 else 0.0
            deltas.append({
                "module": mod, "delta": d, "p_value": float(pval),
                "effect_size": float(es),
                "mean_condition1": float(np.mean(v1)),
                "mean_condition2": float(np.mean(v2)),
            })
        return pd.DataFrame(deltas).set_index("module")

    delta_df = _compare_one(scores_condition1, scores_condition2)
    out = {"delta_modules": delta_df}

    if project_to_eigenspace:
        ref = _load_reference(reference) if (reference is None or isinstance(reference, str)) else reference
        module_order = ref["module_order"]
        n_pcs = ref["n_pcs"]

        loadings = np.zeros((len(module_order), n_pcs))
        for j in range(n_pcs):
            for i, mod in enumerate(module_order):
                loadings[i, j] = ref["loadings"][f"PC{j+1}"].get(mod, 0.0)

        delta_vec = np.array([
            delta_df.loc[mod, "delta"] if mod in delta_df.index else 0.0
            for mod in module_order
        ])
        delta_12d = delta_vec @ loadings
        pc_names = [f"PC{i+1}" for i in range(n_pcs)]
        out["delta_eigenspace"] = pd.Series(delta_12d, index=pc_names, name="delta")
        out["magnitude"] = float(np.linalg.norm(delta_12d))

    top = delta_df.reindex(delta_df["delta"].abs().sort_values(ascending=False).index).head(10)
    out["top_modules"] = top

    delta_vals = delta_df["delta"].values
    n_up = int(np.sum(delta_vals > 0))
    n_down = int(np.sum(delta_vals < 0))
    mean_d = float(np.mean(delta_vals))
    if n_up > 0.7 * len(delta_vals) and mean_d > 0:
        out["disease_type"] = "convergent_activation"
    elif n_down > 0.7 * len(delta_vals) and mean_d < 0:
        out["disease_type"] = "convergent_suppression"
    else:
        out["disease_type"] = "divergent"

    out["infrastructure_regime"] = infrastructure_regime(delta_df)

    if cell_type_column and obs1 is not None and obs2 is not None:
        per_ct = {}
        cell_types = sorted(set(obs1.unique()) & set(obs2.unique()))
        ct_vectors = {}

        for ct in cell_types:
            m1 = obs1 == ct
            m2 = obs2 == ct
            if m1.sum() < 10 or m2.sum() < 10:
                continue
            ct_delta = _compare_one(scores_condition1.loc[m1], scores_condition2.loc[m2])

            if project_to_eigenspace:
                vec = np.array([
                    ct_delta.loc[mod, "delta"] if mod in ct_delta.index else 0.0
                    for mod in module_order
                ])
                ct_12d = vec @ loadings
                ct_vectors[ct] = ct_12d
                eig_vec = ct_12d.tolist()
                mag = float(np.linalg.norm(ct_12d))
            else:
                eig_vec, mag = None, None

            per_ct[ct] = {
                "delta_modules": ct_delta,
                "eigenspace_vector": eig_vec,
                "magnitude": mag,
                "n_condition1": int(m1.sum()),
                "n_condition2": int(m2.sum()),
                "infrastructure_regime": infrastructure_regime(ct_delta),
            }

        if len(ct_vectors) >= 2:
            ct_names = list(ct_vectors.keys())
            cosines = {}
            for i in range(len(ct_names)):
                for j in range(i + 1, len(ct_names)):
                    cosines[f"{ct_names[i]}_vs_{ct_names[j]}"] = cosine(
                        ct_vectors[ct_names[i]], ct_vectors[ct_names[j]]
                    )
            mean_cos = float(np.mean(list(cosines.values())))
            per_ct["_cross_cell_type"] = (
                "convergent" if mean_cos > 0.5 else "anti-convergent" if mean_cos < -0.3 else "divergent"
            )
            per_ct["_cosine_matrix"] = cosines

        out["per_cell_type"] = per_ct

    return out


_IMMUNE = ("monocyte", "macrophage", "t cell", "b cell", "nk", "dendritic", "myeloid",
           "lymphocyte", "immune", "mast")
_STRUCTURAL = ("epithelial", "fibroblast", "endothelial", "hepatocyte", "cardiomyocyte",
               "myocyte", "podocyte", "tubular", "alveolar", "stromal", "pericyte", "adipocyte")


def divergence_score(per_cell_type_results, immune_keywords=_IMMUNE, structural_keywords=_STRUCTURAL):
    """Cosine between mean immune and mean structural eigenspace vectors.

    Negative = divergent, near zero = neutral, positive = convergent.
    """
    immune_vecs, structural_vecs, immune_t, structural_t = [], [], [], []

    for ct_name, data in per_cell_type_results.items():
        if ct_name.startswith("_") or not isinstance(data, dict):
            continue
        v = data.get("eigenspace_vector")
        if v is None:
            continue
        v = np.array(v)
        ct_low = ct_name.lower()
        if any(kw in ct_low for kw in immune_keywords):
            immune_vecs.append(v)
            immune_t.append(ct_name)
        elif any(kw in ct_low for kw in structural_keywords):
            structural_vecs.append(v)
            structural_t.append(ct_name)

    if not immune_vecs or not structural_vecs:
        return {
            "divergence_score": None, "immune_types": immune_t,
            "structural_types": structural_t, "pattern": "insufficient_data",
        }

    im = np.mean(immune_vecs, axis=0)
    st = np.mean(structural_vecs, axis=0)
    score = cosine(im, st)
    pat = "divergent" if score < -0.3 else "convergent" if score > 0.3 else "neutral"
    return {
        "divergence_score": float(score),
        "immune_types": immune_t, "structural_types": structural_t,
        "immune_mean_vector": im.tolist(), "structural_mean_vector": st.tolist(),
        "pattern": pat,
    }
