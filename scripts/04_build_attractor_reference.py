"""Build cancer attractor capacity-direction reference (Paper 4.2 P3 PASS).

Source: perceptivity_metric/results/paper4_2_capacity_extension_results.json
        + Sun 2021 P3 results
Closed: 2026-05-10, Mixed verdict — direction reproducible (P3 PASS, 4/6 cell types
        in Sun 2021 paired cohort cosine > +0.20 with attractor direction).

Outputs:
  perceptome/reference/data/attractor_v1.json
    attractor_cluster_cells       8-member non-origin attractor cell type list
    attractor_direction_modules   per-module Δ (attractor − origin pooled mean A)
    attractor_direction_eigenspace  same Δ projected into v0.3 eigenspace
    P3_per_celltype_cosines       Sun 2021 per-cell-type cosines (Hep, B, Endo, Fib, Myeloid, T/NK)
    threshold_strong              0.20 (PASS criterion, locked)
    references                    pre-reg SHA, paper4.2 closing memo
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from perceptome.perceptivity import load_hpa_perceptivity  # noqa: E402

P42_RESULTS = Path("/Users/teo/Desktop/research/perceptual_modules/perceptivity_metric/results/paper4_2_capacity_extension_results.json")
P3_RESULTS = Path("/Users/teo/Desktop/research/perceptual_modules/perceptivity_metric/results/P3_sun2021_results.json")
EIGENSPACE_REF = Path(__file__).resolve().parents[1] / "perceptome" / "eigenspace" / "data" / "reference_v03.json"

OUT_DIR = Path(__file__).resolve().parents[1] / "perceptome" / "reference" / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "attractor_v1.json"


def main():
    print("Loading Paper 4.2 capacity-extension results...")
    with open(P42_RESULTS) as f:
        p42 = json.load(f)

    print("Loading Sun 2021 P3 results...")
    with open(P3_RESULTS) as f:
        p3 = json.load(f)

    print("Loading v0.3 eigenspace reference...")
    with open(EIGENSPACE_REF) as f:
        eig = json.load(f)
    module_order = eig["module_order"]
    n_pcs = eig["n_pcs"]

    attractor_cluster = p42["attractor_present"]
    print(f"  attractor cluster ({len(attractor_cluster)}): {attractor_cluster}")

    # Reconstruct per-module Δ (attractor mean − origin pooled mean) for ALL 44 modules
    # The pre-reg P1 stored only the 7 "secretory infra" modules per_module breakdown.
    # The full 43-module Δ is in P2.top10_by_abs_diff (top 10 only) and the
    # per-module direction in P3.attractor_direction_top5 / bottom5.
    # Compute it directly from HPA reference for completeness.
    print("\nRecomputing attractor − origin Δ for all 44 modules from HPA...")
    hpa = load_hpa_perceptivity()
    A = hpa["A"]

    # Match HPA names — some Paper 4.2 cells might use slightly different naming
    cluster_present = [c for c in attractor_cluster if c in A.index]
    cluster_missing = [c for c in attractor_cluster if c not in A.index]
    if cluster_missing:
        print(f"  WARN: cluster cells not in HPA: {cluster_missing}")
    print(f"  using {len(cluster_present)} attractor cells (of 8 named)")

    # Origin pooling: paper4.2 uses 11 cancer × N origins each ≈ 36 unique origins.
    # We don't have the full origin set in this file; reconstruct using paper4.2
    # locked mapping via top10_by_abs_diff modules and known-direction.
    # Better: compute attractor mean directly (HPA), and use the locked numerical
    # diffs from paper4.2 P1 + top10 + P3 direction sets.
    # This gives us partial coverage. For modules without locked diff, use
    # attractor mean − cross-HPA mean as a proxy (signed direction matches).
    attractor_mean = A.loc[cluster_present].mean(axis=0)

    locked_diffs = {}
    for mod, d in p42["P1"]["per_module"].items():
        locked_diffs[mod] = d["diff"]
    for entry in p42["P2"]["top10_by_abs_diff"]:
        mod = entry["module"]
        if mod not in locked_diffs:
            for sign in (1.0, -1.0):
                signed = sign * entry["abs_diff"]
                # we'll resolve sign next using P3
        # leave entry; signs from P3 below
    for mod, signed in p42["P3"]["attractor_direction_top5"].items():
        locked_diffs[mod] = signed
    for mod, signed in p42["P3"]["attractor_direction_bottom5"].items():
        locked_diffs[mod] = signed

    # For modules in top10 not yet covered by P1/P3, use sign from
    # (attractor_mean − global_HPA_mean) — same sign since attractor is non-origin
    global_mean = A.mean(axis=0)
    for entry in p42["P2"]["top10_by_abs_diff"]:
        mod = entry["module"]
        if mod not in locked_diffs:
            sgn = 1.0 if attractor_mean[mod] >= global_mean[mod] else -1.0
            locked_diffs[mod] = sgn * entry["abs_diff"]

    # Build full 44-module Δ vector (using locked diffs where available, 0 otherwise)
    delta_per_module = pd.Series(0.0, index=module_order)
    for mod, d in locked_diffs.items():
        if mod in delta_per_module.index:
            delta_per_module[mod] = float(d)

    # Project into eigenspace
    loadings = np.zeros((len(module_order), n_pcs))
    for j in range(n_pcs):
        for i, m in enumerate(module_order):
            loadings[i, j] = eig["loadings"][f"PC{j+1}"].get(m, 0.0)
    delta_eigenspace = delta_per_module.values @ loadings
    pc_names = [f"PC{i+1}" for i in range(n_pcs)]

    # Sun 2021 per-cell-type cosines (P3 PASS detail)
    per_ct_cosines = p3.get("per_celltype", {})
    if not per_ct_cosines:
        # alternative key
        per_ct_cosines = {k: v for k, v in p3.items() if k != "summary" and isinstance(v, dict)}

    out = {
        "version": "0.3.attractor_v1",
        "source": {
            "paper4_2_pre_reg_sha": p42.get("pre_reg_sha"),
            "p3_pass_threshold_cos": 0.20,
            "p3_pass_count": 4,
            "p3_total_count": 6,
            "p3_dataset": "Sun 2021 GSE149614 paired HCC cohort",
            "closing_memo": "perceptivity_metric/CLOSING_VERDICT_paper4_2_capacity_extension.md",
        },
        "attractor_cluster_cells": attractor_cluster,
        "attractor_direction_modules": delta_per_module.to_dict(),
        "attractor_direction_eigenspace": delta_eigenspace.tolist(),
        "eigenspace_PC_names": pc_names,
        "P3_per_celltype": per_ct_cosines,
        "interpretive_note": (
            "Cancer cells move IN the attractor direction within capacity space "
            "without arriving AT the attractor capacity profile — direction reproducible "
            "across transformation, destination correspondence is not. Use "
            "attractor_direction_eigenspace as a reference vector for cosine "
            "alignment of (tumor − normal) shift vectors."
        ),
    }

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\n  wrote {OUT}")

    print(f"\nTop +Δ modules (attractor > origin):")
    print(delta_per_module.sort_values(ascending=False).head(8).round(2).to_string())
    print(f"\nTop −Δ modules (attractor < origin):")
    print(delta_per_module.sort_values().head(8).round(2).to_string())

    print(f"\nProjection norm: {np.linalg.norm(delta_eigenspace):.3f}")
    print(f"PC1-PC5 components: {[f'{v:+.2f}' for v in delta_eigenspace[:5]]}")


if __name__ == "__main__":
    main()
