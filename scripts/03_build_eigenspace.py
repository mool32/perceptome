"""Build the v0.3 reference eigenspace (12-PC from 154 × 44 HPA).

Method matches Paper 3 v0.2 procedure (uses READINESS, R = core_genes):
  1. Load 154 × 44 R matrix from HPA perceptivity reference (= mean log1p(nCPM)
     of core_genes per cell type per module). This matches the v0.2 paper3
     pipeline where score_modules(gene_set='core', method='mean_raw') was the default.
  2. Z-score each module across cell types.
  3. Correlation matrix → eigendecomposition.
  4. Kaiser criterion (eigenvalue > 1).
  5. Bootstrap 100 resamples for per-PC stability.
  6. Save to perceptome/eigenspace/data/reference_v03.json.

A note on R vs A as eigenspace base:
  R-based eigenspace (this file, matches v0.2): "what machinery is present"
                                                projection works for any gene_set='core' scores.
  A-based eigenspace (alternative): "what is currently engaged" — would require
                                    projecting gene_set='activity' scores.
  For v0.3 default we keep R-based for paper3 compatibility. An A-based eigenspace
  is reserved for future v0.x as `pct.eigenspace.activity_reference()` if the
  user community asks for it.
"""

import sys
from pathlib import Path
import json

import numpy as np
import pandas as pd

# Make the in-tree package importable
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from perceptome.perceptivity import load_hpa_perceptivity  # noqa: E402
from perceptome.eigenspace.rebuild import rebuild  # noqa: E402

OUT_DIR = Path(__file__).resolve().parents[1] / "perceptome" / "eigenspace" / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "reference_v03.json"


def main():
    print("Loading HPA perceptivity reference (154 × 44)...")
    ref = load_hpa_perceptivity()
    R = ref["R"]
    print(f"  R matrix shape (readiness — core_genes): {R.shape}")
    print(f"  NaN count: {R.isna().sum().sum()} (should be ~0)")

    if R.isna().any().any():
        print("  filling NaN with 0 (cells with no core genes detected)")
        R = R.fillna(0.0)

    print("\nRebuilding eigenspace from R (Kaiser cutoff, 100 bootstraps)...")
    out = rebuild(R, kaiser=True, bootstrap_n=100, seed=42)

    print(f"\n  n_kaiser (eigenvalues > 1): {out['n_kaiser']}")
    print(f"  Top 5 eigenvalues: {[f'{v:.2f}' for v in out['eigenvalues'][:5]]}")
    print(f"  Bottom 3 eigenvalues: {[f'{v:.4f}' for v in out['eigenvalues'][-3:]]}")
    print(f"  Variance explained PC1-PC5: {[f'{v:.3f}' for v in out['var_explained'][:5]]}")
    print(f"  Cumulative var PC1-PC{out['n_kaiser']}: {sum(out['var_explained'][:out['n_kaiser']]):.3f}")
    print(f"  Spectral fit α = {out['spectral_fit']['alpha']:.3f}, R² = {out['spectral_fit']['R2']:.3f}")

    print(f"\n  Bootstrap stability (per PC):")
    for pc, s in sorted(out["stability"].items(), key=lambda kv: int(kv[0][2:])):
        flag = "stable" if s > 0.6 else "probable" if s > 0.4 else "exploratory"
        print(f"    {pc}: {s:.3f}  ({flag})")

    out["version"] = "0.3"
    out["catalog_version"] = "0.3"
    out["construction_method"] = "z-score across cell types, correlation eigendecomp, Kaiser cutoff"

    pc_confidence = {}
    for pc, s in out["stability"].items():
        pc_confidence[pc] = "stable" if s > 0.6 else "probable" if s > 0.4 else "exploratory"
    out["pc_confidence"] = pc_confidence

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\n  wrote {OUT} ({OUT.stat().st_size / 1024:.1f} KB)")

    # Compare with v0.2 (43-module) eigenspace if available
    v02_path = Path("/Users/teo/Desktop/research/perceptual_modules/tool/perceptome/perceptome/data/reference_eigenspace.json")
    if v02_path.exists():
        with open(v02_path) as f:
            v02 = json.load(f)
        print(f"\nComparison with v0.2 (43-module, paper3 release):")
        print(f"  v0.2 n_pcs: {v02.get('n_pcs')}, modules: {v02.get('n_modules')}")
        print(f"  v0.3 n_pcs: {out['n_kaiser']}, modules: {out['n_modules']}")
        print(f"  Top 5 eigenvalues v0.2: {[f'{v:.2f}' for v in v02.get('eigenvalues', [])[:5]]}")

    print("\nValidation: project an HPA cell type back into the new eigenspace")
    from perceptome.eigenspace.project import project, _load_reference
    _load_reference.cache_clear()
    sample_scores = R.iloc[:3]  # first 3 cell types
    coords = project(sample_scores)
    print(f"  Sample projections (PC1-PC3):")
    print(coords["coordinates"][["PC1", "PC2", "PC3"]].round(3).to_string())


if __name__ == "__main__":
    main()
