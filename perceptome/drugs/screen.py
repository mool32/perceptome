"""Activity-layer screen — the validated drug operation (Block 5 v1.2 procedure).

For a user perturbation dataset (signatures × genes), test whether a class of
drugs (or a single drug) engages a specific module's TF-target panel above
what a background pool of unrelated drugs would predict by chance.

This is the operation Paper 4.1 §1 validated. Three other obvious-looking
operations are FALSIFIED — see PAPER4_1_CANONICAL_RESULTS.md, do not use:
  • cosine matching of drug positions to disease vectors in eigenspace
  • mechanism deconvolution from drug centroid geometry
  • readiness-layer screening (use activity layer; readiness conflates TF + targets)
"""

import numpy as np
import pandas as pd

from ..catalog import get_genes


def activity_layer_screen(
    adata,
    pert_col,
    test_perturbations,
    panels="all_validated",
    background="all_other",
    n_perm=10000,
    seed=42,
    expression_layer=None,
):
    """Class-level activity-layer test against background null.

    Parameters
    ----------
    adata : AnnData
        Signatures (or cells) × genes. .X expected to be drug perturbation
        signatures (e.g., LINCS L1000 z-scores) or log-normalized expression.
    pert_col : str
        adata.obs column naming the perturbation (drug) for each row.
    test_perturbations : str | list | dict
        - str       : single drug name
        - list[str] : list of drug names treated as one class
        - dict      : {class_name: [drug_names]} for multi-class panel
    panels : 'all_validated' | 'rescues_only' | 'controls_only' | list of (class, module, sign)
        Which (class, module, sign) triples to test. 'all_validated' = the 9
        anchors. Custom triples allowed for exploratory work but lose Paper 4.1
        validation guarantees.
    background : 'all_other' | list[str]
        How to define the null pool. 'all_other' = all perturbations in adata
        not in test_perturbations. Or pass an explicit list of drug names.
    n_perm : int
        Permutation count for null distribution.
    seed : int
        RNG seed.
    expression_layer : str | None
        adata.layers key to score on. None uses adata.X.

    Returns
    -------
    DataFrame with one row per (test_class, panel_class, panel_module) tested:
        observed_z          mean activity score across class members
        null_median, null_mean, null_std
        p_one_sided
        q_BH                Benjamini–Hochberg FDR across all rows
        verdict             'PASS' (q<0.10) | 'FAIL' | 'NA'
        n_class_members
        expected_sign

    Notes
    -----
    The procedure exactly mirrors Block 5 v1.2 (step28): one-sided permutation
    against background null, BH FDR. PASS at q<0.10 matches the locked Paper 4.1
    pre-reg threshold.
    """
    from .anchors import ANCHORS

    rng = np.random.default_rng(seed)

    if isinstance(test_perturbations, str):
        test_classes = {test_perturbations: [test_perturbations]}
    elif isinstance(test_perturbations, list):
        test_classes = {"_user_class": test_perturbations}
    elif isinstance(test_perturbations, dict):
        test_classes = test_perturbations
    else:
        raise TypeError("test_perturbations must be str | list | dict")

    if isinstance(panels, str):
        if panels == "all_validated":
            triples_df = ANCHORS
        elif panels == "rescues_only":
            triples_df = ANCHORS[ANCHORS["role"] == "validated_rescue"]
        elif panels == "controls_only":
            triples_df = ANCHORS[ANCHORS["role"] == "positive_control"]
        else:
            raise ValueError(f"panels string must be one of all_validated|rescues_only|controls_only, got {panels!r}")
        triples = list(zip(triples_df["class"], triples_df["module"], triples_df["expected_sign"]))
    else:
        triples = list(panels)

    pert = adata.obs[pert_col].astype(str)
    test_drugs_all = set()
    for drugs in test_classes.values():
        test_drugs_all.update(drugs)

    if isinstance(background, str) and background == "all_other":
        bg_drugs = sorted(set(pert.unique()) - test_drugs_all)
    elif isinstance(background, list):
        bg_drugs = list(background)
    else:
        raise TypeError("background must be 'all_other' or a list of drug names")

    if len(bg_drugs) < 30:
        raise ValueError(
            f"Background pool too small ({len(bg_drugs)}). Need ≥30; ideally ~1000 "
            "as in the validated Paper 4.1 procedure."
        )

    # Score each row × module (mean activity-layer score across found genes)
    module_names_needed = sorted({m for _, m, _ in triples})
    module_to_genes = {m: [g for g in get_genes(m, "activity") if g in adata.var_names]
                       for m in module_names_needed}

    sig_module_scores = {}
    for m in module_names_needed:
        genes = module_to_genes[m]
        if not genes:
            sig_module_scores[m] = pd.Series(np.nan, index=adata.obs_names)
            continue
        if expression_layer:
            X_sub = adata[:, genes].layers[expression_layer]
        else:
            X_sub = adata[:, genes].X
        if hasattr(X_sub, "toarray"):
            X_sub = X_sub.toarray()
        sig_module_scores[m] = pd.Series(
            np.mean(np.asarray(X_sub, dtype=float), axis=1),
            index=adata.obs_names,
        )

    # Per-drug × module aggregate (mean across this drug's signatures)
    drug_module_scores = {}
    for drug, mask in pert.groupby(pert).groups.items():
        drug_module_scores[drug] = {
            m: float(np.nanmean(sig_module_scores[m].loc[list(mask)]))
            for m in module_names_needed
        }

    # Background drug × module matrix
    bg_module_matrix = np.array([
        [drug_module_scores.get(d, {}).get(m, np.nan) for m in module_names_needed]
        for d in bg_drugs
    ])
    mod_to_col = {m: i for i, m in enumerate(module_names_needed)}

    rows = []
    for test_class, drugs in test_classes.items():
        class_present = [d for d in drugs if d in drug_module_scores]
        n_class = len(class_present)
        if n_class == 0:
            for cls, mod, sign in triples:
                rows.append({
                    "test_class": test_class,
                    "panel_class": cls,
                    "panel_module": mod,
                    "expected_sign": sign,
                    "observed_z": np.nan,
                    "null_median": np.nan,
                    "null_mean": np.nan,
                    "null_std": np.nan,
                    "p_one_sided": np.nan,
                    "n_class_members": 0,
                    "verdict": "NA_no_drugs_found",
                })
            continue

        for cls, mod, sign in triples:
            obs_vals = [drug_module_scores[d][mod] for d in class_present
                        if not np.isnan(drug_module_scores[d][mod])]
            if not obs_vals:
                rows.append({
                    "test_class": test_class, "panel_class": cls, "panel_module": mod,
                    "expected_sign": sign, "observed_z": np.nan,
                    "null_median": np.nan, "null_mean": np.nan, "null_std": np.nan,
                    "p_one_sided": np.nan, "n_class_members": n_class,
                    "verdict": "NA_module_genes_missing",
                })
                continue
            observed = float(np.mean(obs_vals))
            bg_col = bg_module_matrix[:, mod_to_col[mod]]
            bg_clean = bg_col[~np.isnan(bg_col)]
            if len(bg_clean) < 30:
                rows.append({
                    "test_class": test_class, "panel_class": cls, "panel_module": mod,
                    "expected_sign": sign, "observed_z": observed,
                    "null_median": np.nan, "null_mean": np.nan, "null_std": np.nan,
                    "p_one_sided": np.nan, "n_class_members": n_class,
                    "verdict": "NA_bg_pool_too_small",
                })
                continue

            idx = rng.integers(0, len(bg_clean), size=(n_perm, n_class))
            null_means = bg_clean[idx].mean(axis=1)
            if sign > 0:
                p = float((null_means >= observed).sum() + 1) / (n_perm + 1)
            else:
                p = float((null_means <= observed).sum() + 1) / (n_perm + 1)

            rows.append({
                "test_class": test_class, "panel_class": cls, "panel_module": mod,
                "expected_sign": sign,
                "observed_z": observed,
                "null_median": float(np.median(null_means)),
                "null_mean": float(null_means.mean()),
                "null_std": float(null_means.std()),
                "p_one_sided": p,
                "n_class_members": n_class,
            })

    df = pd.DataFrame(rows)

    # BH FDR across non-NA rows
    valid = df["p_one_sided"].notna()
    pv = df.loc[valid, "p_one_sided"].values
    if len(pv) > 0:
        order = np.argsort(pv)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, len(pv) + 1)
        q = pv * len(pv) / ranks
        # monotonicity (BH)
        q_sorted = q[order]
        q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
        q_final = np.empty_like(q)
        q_final[order] = q_sorted
        df.loc[valid, "q_BH"] = np.minimum(q_final, 1.0)
    else:
        df["q_BH"] = np.nan

    def _verdict(row):
        if pd.isna(row.get("p_one_sided")):
            return row.get("verdict", "NA")
        return "PASS" if row["q_BH"] < 0.10 else "FAIL"
    df["verdict"] = df.apply(_verdict, axis=1)

    return df
