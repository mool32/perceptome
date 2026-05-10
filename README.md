# perceptome

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python: ≥3.8](https://img.shields.io/badge/Python-≥3.8-blue.svg)](https://www.python.org)
[![Tests: 42/42](https://img.shields.io/badge/tests-42%2F42%20passing-green.svg)](#tests)

**Cellular perception analysis toolkit.** 44 transcriptional signaling modules organized into a 9-PC eigenspace derived from 154 Human Protein Atlas cell types. Score any RNA-seq dataset, project into the canonical perceptual coordinate system, and answer three categories of question:

1. **Geometry** — *where does this cell live in module space?*
2. **Perceptivity** — *what can this cell do?* (capacity, headroom, saturation)
3. **Validity** — *did this measurement reach biology, or hit a confound?*

> **Companion papers:** Paper 3 (perceptome eigenspace) — accompanies the framework. Paper 4.1 (drugs, 6 validated rescues + 11 falsifications) — defines drug-pharmacology scope. Paper 4.2 (cancer attractor capacity-direction). Substrate-series Papers 4.3-4.8 (immune, muscle, epithelial, organoid) define operational scope of capacity-floor predictor.

**New users start here:**

- 📓 [`examples/tutorial.ipynb`](examples/tutorial.ipynb) — executable end-to-end tutorial on real public data (PBMC3K). Runs in <2 min, ships with output cells populated.
- 📖 [`docs/QUICKSTART.md`](docs/QUICKSTART.md) — written form of the tutorial with extra notes on common pitfalls
- 🎯 [`docs/SCOPE.md`](docs/SCOPE.md) — what the framework answers, what it doesn't (8 falsified operations explicitly listed)
- 🧮 [`docs/EIGENSPACE_EVOLUTION.md`](docs/EIGENSPACE_EVOLUTION.md) — relationship between v0.1 12-PC and v0.2 9-PC eigenspaces (read if comparing tool output to Paper 3)

---

## Installation

```bash
pip install git+https://github.com/mool32/perceptome.git
# or, for development:
git clone https://github.com/mool32/perceptome.git
cd perceptome
pip install -e .
```

Requires Python ≥ 3.8 with scanpy, anndata, numpy, pandas, scipy, scikit-learn.

## Quickstart

```python
import scanpy as sc
import perceptome as pct

adata = sc.read_h5ad("your_data.h5ad")  # log-normalized

# 1. GEOMETRY — score 44 modules per cell, project into 9-PC eigenspace
scores = pct.score_modules(adata, method="mean_raw")["scores"]
coords = pct.project(scores)["coordinates"]
refs = pct.compare_to_references(coords)  # cosine vs cancer attractor

# 2. PERCEPTIVITY — what can each cell type do?
R = pct.score_readiness(adata)            # mean log expression of core_genes
A = pct.score_activity(adata)             # mean log expression of TF-target genes
perc = pct.compute_perceptivity(
    R, A,
    cell_type=adata.obs["cell_type"],
    cell_class=adata.obs["tissue"],
)
# per-cell-type DataFrame with R, A, C, headroom, BS, GS, I, spec_quadrant

# Predict architecture engagement a priori for a known operation
pred = pct.predict_engagement(
    starting_cell_type="cardiomyocytes",
    operation_modules=["HSF1", "UPR-ATF6", "mTOR"],
)
# → capacity_floor (saturated_blocked_up | capacious | intermediate)
#   predicted_direction (up_blocked | up_possible | …)

# 3. VALIDITY — did the perturbation reach biology, or hit a confound?
sc_obj = pct.validate_perturbation(
    adata, condition_col="treatment",
    perturbation_value="drug", control_value="ctrl",
)
print(pct.scorecard(sc_obj))
# Overall: PASS | MIXED | ARTIFACT | INCONCLUSIVE
# Checks:  ✓ random_200    ✓ housekeeping    ✓ cell_cycle
```

## What's in the box

| Layer | Function | What it answers |
|---|---|---|
| **catalog** | `pct.list_modules()`, `pct.get_genes(module, gene_set)` | What modules and gene panels are defined? |
| **score** | `pct.score_modules`, `pct.score_readiness`, `pct.score_activity` | Module activity per cell from expression |
| **eigenspace** | `pct.project`, `pct.rebuild` | Project into 9-PC HPA-derived perceptual space |
| **perceptivity** | `pct.compute_perceptivity`, `pct.predict_engagement`, `pct.capacity_floor`, `pct.load_hpa_perceptivity` | 5-vector capacity / engagement / specialization per cell type |
| **compare** | `pct.compare_conditions`, `pct.compare_to_references`, `pct.infrastructure_regime`, `pct.divergence_score` | Disease vs control deltas, regime classification |
| **reference** | `pct.load_attractor_direction`, `pct.attractor_alignment` | Cancer capacity-direction (Paper 4.2 P3) |
| **validity** | `pct.validate_perturbation`, `pct.scorecard` | Random-200 / housekeeping / cell-cycle nulls |
| **drugs** | `pct.drug_anchors`, `pct.activity_layer_screen` | 9-anchor narrow validated drug operation |
| **network** | `pct.compute_network`, `pct.module_heterogeneity` | Module co-variation + rich-club analysis |

## Reference data bundled

- **44-module catalog** with core / activity / sensor / cascade / TF / feedback gene lists per module (perceptome v0.3 catalog: v0.1 + NPAS4)
- **154 × 44 HPA perceptivity reference** (R, A, C, headroom matrices, precomputed)
- **9-PC HPA eigenspace** (12-PC v0.2 superseded by 9-PC v0.3 from 154 × 44 R matrix; PC1-PC6 stable per bootstrap)
- **Cancer attractor capacity-direction** (Paper 4.2 P3 PASS, 8-cell-type cluster, locked Δ vector)
- **9-anchor drug reference table** (Paper 4.1 validated mechanism-pathway pairs, with full validation chain — Block 5 v1.2 z + q + Block 2 holdout drugs + Block 4 LD)
- **Validity panels**: housekeeping (Eisenberg & Levanon 2013), cell cycle (Tirosh 2016)

**Removed in v0.2:** the 19,804-compound CMap drug eigenspace projection from v0.1 was a pre-computed asset for drug-disease cosine matching, an operation that Paper 4.1 falsified in three independent formulations. Replaced by the 9-anchor `drug_anchors` table + `activity_layer_screen` operation, which is what Paper 4.1 actually validated.

**Pinned for backward compatibility:** v0.2 disease (RA/AD/IPF/DKD) and aging (inflammaging/collapse) reference vectors are retained as `*_v02pinned.json` (12-PC v0.2 space, incompatible with v0.3 9-PC eigenspace). Recompute on v0.3 scheduled for v0.2.1.

## Methodology

- **`mean_raw`** is the default scoring method; safe across datasets. **`scanpy_corrected`** (background-corrected) recommended when comparing cells with very different transcriptome complexity (e.g., proliferating vs post-mitotic). Switch was triggered by Paper 4.5 v1.2 amendment after random-200 ARTIFACT FAIL detection in adult-intestinal-epithelium pipeline.
- **Module activity = readiness** (R, core_genes) and **engagement** (A, TF-target genes) are different things. R-A correlation across 43 modules is ρ ≈ 0.75; CREB and NPAS4 are systematic exceptions (post-translational regulation, rapid protein turnover). Use both for any module marked `dissociation_risk: HIGH` in the catalog.
- **The 9-PC v0.3 eigenspace is fixed reference**; new datasets project into it without refitting. Paper 3 v0.2 used 12-PC, 43-module — superseded by v0.3 (44 modules with NPAS4, recomputed eigenspace, PC1-PC6 stable vs PC1-PC3 in v0.2).

### The capacity-floor predictor (closed paper4.5 + 4.6 + 4.7 + 4.8)

For any cell type × module pair from HPA:

- `A_baseline > 4.5` → **saturated**, cannot ramp UP further (downward suppression by specific signaling, e.g. retinoic acid → UPR, IS possible — predictor is upward-asymmetric)
- `A_baseline < 2.5` → **capacious**, ramp possible; magnitude is operation-determined
- `2.5 ≤ A ≤ 4.5` → intermediate

This is **factor 1** of the two-factor framework (`engagement = capacity × operation_intensity`). Factor 2 (operation intensity) requires context-specific modeling outside the tool.

### Drug perturbation analysis — narrow scope

Paper 4.1 closed 2026-05-09 with **6 surviving validated findings + 11 pre-registered falsifications**. The validated drug operation is **activity-layer scoring of TF-target panels** for 9 specific (class, module) anchors — see `pct.drug_anchors()` for the table with full validation chain.

**Falsified drug operations — `perceptome` does not support these and you should not use them:**
- Drug-disease cosine matching in eigenspace (3 formulations)
- Drug-class mechanism deconvolution from panel geometry
- Readiness-layer timescale rescue (24h vs 6h)
- TF-autoregulation as 1st-class layer
- Equilibrium-instrument linear-response framing
- Clean perceptome signature predicting FDA approval beyond selectivity
- Snapshot perceptome blind to pulsatile dynamics

See `PAPER4_1_CANONICAL_RESULTS.md` and `killed_hypotheses.csv` in the source tree for the full audit trail.

## Tests

```bash
pip install pytest
pytest tests/ -v
# 42/42 passing in <1 second
```

Coverage: catalog, scoring, perceptivity (capacity floor + a-priori predictions reproducing paper4.5 / 4.4 results), eigenspace projection, attractor reference, validity scorecard, drug anchors + activity_layer_screen.

## Citation

```bibtex
@software{perceptome2026,
  author = {Spiro, Theodor},
  title = {perceptome: Cellular Perception Analysis Toolkit},
  version = {0.2.0},
  year = {2026},
  url = {https://github.com/mool32/perceptome}
}
```

## Versioning

This is **v0.2.0**. Major changes vs v0.1.0:
- 44-module catalog (NPAS4 added per Paper 4 Table S4)
- 9-PC eigenspace recomputed from 154 × 44 R matrix
- New `pct.perceptivity` subpackage (capacity-floor predictor + two-factor framework)
- New `pct.validity` subpackage (random-200/housekeeping/cell-cycle scorecard)
- New `pct.reference.attractor` (Paper 4.2 P3 capacity-direction)
- New `pct.drugs` subpackage (9-anchor table + activity_layer_screen) replacing 19K-compound vectors
- Restructured into per-concern subpackages (was single flat module)

See [`CHANGELOG.md`](CHANGELOG.md) for full diff and [`ROADMAP.md`](ROADMAP.md) for v0.3+ direction.

## License

MIT — see [`LICENSE`](LICENSE).

## Author

Theodor Spiro · ORCID [0009-0004-5382-9346](https://orcid.org/0009-0004-5382-9346) · tspiro@vaika.org · Vaika Inc., East Aurora, NY, USA
