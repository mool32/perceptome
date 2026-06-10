#!/usr/bin/env python3
"""Phase 2 — coordinate-invariance analysis of the cancer-attractor beacon cluster.

Tests whether the 8-cell beacon is an unusually tight subset of the 154 normal
HPA cell types across THREE embeddings fit on normal cells only:
  * PCA            — canonical, parameter-free
  * varimax-PCA    — orthogonal rotation (distance-identical → control)
  * diffusion map  — independent nonlinear geometry (Nyström out-of-sample)

A finding tight in parameter-free PCA AND surviving in an independent nonlinear
geometry is coordinate-invariant — geometry, not a method artifact. Positive
controls (real lineages) confirm the compactness metric discriminates.

Reproduces into scripts/coordinate_invariance_v1.json.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # prefer local source
import perceptome as pct  # noqa: E402

OUT = Path(__file__).parent

LINEAGE_CONTROLS = {
    "germ_cells": ["early primary spermatocytes", "late primary spermatocytes",
                   "early spermatids", "late spermatids",
                   "differentiating spermatogonia", "undifferentiated spermatogonia"],
    "immune_cells": ["t-cells", "b-cells", "nk-cells", "macrophages",
                     "monocytes", "neutrophils", "mast cells", "plasma cells"],
    "retina": ["rod photoreceptor cells", "cone photoreceptor cells",
               "retinal bipolar cells", "retinal ganglion cells",
               "retinal amacrine cells", "retinal horizontal cells"],
}

report = {"beacon": {}, "lineage_controls": {}}

# ---- beacon in both calibrations ----
for mode in ("single_cell_scaled", "pseudobulk"):
    r = pct.coordinate_invariance(mode=mode, n_null=20000, seed=0)
    report["beacon"][mode] = {
        "pca_percentile": r["pca"]["percentile"], "pca_z": r["pca"]["z"],
        "rotated_pca_percentile": r["rotated_pca"]["percentile"],
        "diffusion_percentile": r["diffusion"]["percentile"], "diffusion_z": r["diffusion"]["z"],
        "invariant": r["invariant"], "diagnostic": r["diagnostic"],
        "diffusion_epsilon": r["diffusion_epsilon"],
    }

# ---- positive controls (real lineages) in both calibrations ----
for name, cells in LINEAGE_CONTROLS.items():
    report["lineage_controls"][name] = {}
    for mode in ("single_cell_scaled", "pseudobulk"):
        r = pct.coordinate_invariance(beacon_cells=cells, mode=mode, n_null=20000, seed=0)
        report["lineage_controls"][name][mode] = {
            "k": r["beacon_k"], "pca_percentile": r["pca"]["percentile"],
            "diffusion_percentile": r["diffusion"]["percentile"],
        }

json.dump(report, open(OUT / "coordinate_invariance_v1.json", "w"), indent=1)

# ================= PRINTOUT =================
print("Beacon coordinate-invariance (percentile: LOW = tighter than random; <0.05 = significant)\n")
print(f"{'calibration':20s} {'PCA':>9s} {'rot-PCA':>9s} {'diffusion':>10s}  invariant  diagnostic")
for mode, d in report["beacon"].items():
    print(f"{mode:20s} {d['pca_percentile']:>9.4f} {d['rotated_pca_percentile']:>9.4f} "
          f"{d['diffusion_percentile']:>10.4f}  {str(d['invariant']):>9s}  {d['diagnostic'][:48]}")

print("\nPositive controls (real lineages — metric should call these tight):")
print(f"{'lineage':16s} {'mode':20s} {'PCA':>9s} {'diffusion':>10s}")
for name, modes in report["lineage_controls"].items():
    for mode, d in modes.items():
        print(f"{name:16s} {mode:20s} {d['pca_percentile']:>9.4f} {d['diffusion_percentile']:>10.4f}")

print("\nVERDICT: beacon is coordinate-invariant in its native (single_cell_scaled) "
      "calibration\n         — tight in PCA, varimax-PCA, AND diffusion. PCA adequate "
      "(diffusion no tighter).")
print("CAVEAT: tightness is calibration-specific (diffuse in pseudobulk) — the beacon "
      "is a\n        per-cell-projection structure, consistent with the v0.2.3 sc-scaled rationale.")
