# Quickstart — your first analysis with perceptome

You have an scRNA-seq dataset and a question about cellular signaling. This page walks you through the four most common workflows in under 100 lines of code.

## Install

```bash
pip install git+https://github.com/mool32/perceptome.git
# or, for development:
git clone https://github.com/mool32/perceptome.git
cd perceptome
pip install -e .
```

Requires Python ≥3.8 with scanpy, anndata, numpy, pandas, scipy, scikit-learn.

## Setup once

```python
import scanpy as sc
import perceptome as pct
import pandas as pd

adata = sc.read_h5ad("your_data.h5ad")  # log-normalized expression in .X
# Make sure adata.obs has at least one categorical column for cell types,
# and a condition / treatment column if you want to do perturbation analysis.
```

The catalog has **44 modules** — see `pct.list_modules()`. Use `pct.get_genes('NF-κB', 'core')` to inspect any module's gene panel.

---

## Workflow 1 — Geometry: where do my cells live in module space?

The 9-PC HPA-derived eigenspace is the canonical perceptual coordinate system. PC1 is preserved from Paper 3 (perception breadth, myeloid+ to germ−); PC4 carries the cancer convergence direction (Paper 4.2).

```python
scores = pct.score_modules(adata, method="mean_raw")["scores"]
# scores: DataFrame (cells × 44 modules)

coords = pct.project(scores)["coordinates"]
# coords: DataFrame (cells × 9 PCs)

# Add coords back to your AnnData if you want to plot them with scanpy
import numpy as np
adata.obsm["X_perceptome"] = coords.values
sc.pl.scatter(adata, basis="perceptome", color="cell_type")
```

**When to use which scoring method:**
- `mean_raw` (default) — safe for cross-dataset comparisons. Use this unless you have a specific reason not to.
- `scanpy_corrected` — switch to this when comparing cells with very different transcriptome complexity (proliferating vs post-mitotic, healthy vs heavily damaged). Catches the baseline-shift artifact that broke an early Paper 4.5 epithelial analysis.
- `mean_zscore` — within-dataset only. **Never** compare two independently z-scored datasets — it destroys cross-dataset signal.

---

## Workflow 2 — Perceptivity: what can my cells do?

The capacity layer answers a different question from geometry. Geometry tells you *where* a cell sits in module space. Perceptivity tells you *what it can do from that position* — which modules have ramp room (capacity) and which are saturated.

```python
# Score readiness (R) and activity (A) separately
R = pct.score_readiness(adata)   # mean log expression of core_genes per module
A = pct.score_activity(adata)    # mean log expression of TF-target genes per module

perc = pct.compute_perceptivity(
    R, A,
    cell_type=adata.obs["cell_type"],
    cell_class=adata.obs["tissue"],   # for within-class specialization (BS)
)
# perc["per_cell_type"] — DataFrame with R, A, C, BS, GS, I, spec_quadrant per cell type
# perc["per_module"] — dict of cell-type × module matrices for R, A, C, headroom, capacity_floor
```

### A-priori prediction without your data

The HPA reference is bundled. You can ask the capacity-floor predictor a question before running an experiment:

```python
# Will hepatocytes ramp UPR-ATF6 under partial-hepatectomy regeneration?
pred = pct.predict_engagement(
    starting_cell_type="hepatocytes",
    operation_modules=["UPR-ATF6", "HSF1", "mTOR"],
)
print(pred[["A_baseline", "headroom", "capacity_floor", "predicted_direction"]])

#           A_baseline  headroom        capacity_floor       predicted_direction
# UPR-ATF6     6.32      0.96     saturated_blocked_up  up_blocked_down_possible
# HSF1         3.22      1.68             intermediate    intermediate_uncertain
# mTOR         4.88      0.64     saturated_blocked_up  up_blocked_down_possible
```

This is **factor 1** of the two-factor framework (`engagement = capacity × operation_intensity`). Factor 2 (operation intensity) requires context-specific modeling outside the tool. The capacity-floor predictor is **upward-asymmetric**: A_baseline > 4.5 means cells cannot ramp UP further, but specific signaling perturbations (e.g. retinoic acid → UPR) CAN suppress saturated modules. Closed across paper4.5 + 4.6 + 4.7 + 4.8 (2026-05-10).

---

## Workflow 3 — Perturbation analysis: drug vs control

You ran an experiment with a perturbation (drug, knockout, stimulation) and a control. **Run the validity scorecard first** — it catches the most common analysis-breaking artifacts before you trust any module-level effect.

```python
sc_obj = pct.validate_perturbation(
    adata,
    condition_col="treatment",
    perturbation_value="drug",
    control_value="ctrl",
)
print(pct.scorecard(sc_obj))

# Overall: PASS
# Checks:
#   ✓ random_200      value=-0.0042  (< 0.10)  → PASS
#       ↳ baseline-shift artifact detector
#   ✓ housekeeping    value=+0.0123  (< 0.20)  → PASS
#       ↳ technical normalization sanity
#   ✓ cell_cycle      value=+0.4521  (> 0.30)  → PASS
#       ↳ positive control; FAIL ⇒ perturbation may not reach biology
```

| Overall | What it means |
|---|---|
| `PASS` | All three nulls pass — interpret module effects normally |
| `MIXED` | Some nulls fail — interpret with caution |
| `ARTIFACT` | random_200 fails — baseline shift confounds the analysis. Retry with `score_method="scanpy_corrected"` (paper4.5 v1.2 fix) |
| `INCONCLUSIVE` | cell_cycle fails — perturbation may not have reached biology, module signals may be noise |

Once validity is OK, compute the per-module delta + eigenspace shift:

```python
drug_scores = pct.score_modules(adata[adata.obs.treatment == "drug"])["scores"]
ctrl_scores = pct.score_modules(adata[adata.obs.treatment == "ctrl"])["scores"]

delta = pct.compare_conditions(
    drug_scores, ctrl_scores,
    project_to_eigenspace=True,
)
print(delta["top_modules"])           # 10 most-changed modules
print(delta["delta_eigenspace"])      # 9-PC perturbation vector
print(delta["infrastructure_regime"]) # supply_chain | firefighting | collapse | unsupported
```

### Cancer attractor alignment (Paper 4.2)

If your perturbation is a transformation (e.g., tumor vs adjacent normal in paired cohort), test whether the shift aligns with the cancer-transformation direction:

```python
shift_vector = delta["delta_modules"]["delta"]   # 44 modules
result = pct.attractor_alignment(shift_vector, mode="modules")
print(result)
# {'cosine': 0.41, 'passes_p3_threshold': True, 'threshold': 0.20}
# Sun 2021 paired HCC PASSed at 4/6 cell types > +0.20.
```

---

## Workflow 4 — Drug perturbation analysis (narrow validated scope)

**Read this before using.** Paper 4.1 closed 2026-05-09 with **6 surviving validated drug findings + 11 pre-registered falsifications**. The tool supports the validated operation only — see `docs/SCOPE.md` for the full audit trail.

What works: activity-layer scoring of TF-target panels for 9 specific (class, module) anchors, against a background null.

What does NOT work and the tool will not help you do it:
- Drug-disease cosine matching in eigenspace (3 formulations falsified)
- Drug-class mechanism deconvolution from panel geometry
- Readiness-layer timescale rescue
- TF-autoregulation as 1st-class layer
- Equilibrium-instrument linear-response framings
- Clean-signature → FDA approval prediction
- Snapshot perceptome blind to pulsatile dynamics

### View the 9 validated anchors

```python
anchors = pct.drug_anchors()
print(anchors[["class", "module", "expected_sign", "role", "block5_z", "block5_q"]])

#         class       module  expected_sign             role  block5_z  block5_q
#          MEKi     ERK/MAPK             -1 validated_rescue   -1.026    0.0002
#   Proteasomei     UPR-PERK              1 validated_rescue   +1.209    0.0002
#          CDKi   Cell Cycle             -1 validated_rescue      NaN       NaN
#         EGFRi     ERK/MAPK             -1 validated_rescue      NaN       NaN
#         PI3Ki    PI3K/PTEN              1 validated_rescue      NaN       NaN
#          IKKi        NF-κB             -1 validated_rescue      NaN       NaN
#        HSP90i         HSF1              1 positive_control   +2.149    0.0002
#        HSP90i     UPR-ATF6              1 positive_control   +0.736    0.0004
#        HSP90i         NRF2              1 positive_control   +0.510    0.038
```

### Run the screen on your data

```python
result = pct.activity_layer_screen(
    adata,
    pert_col="drug_name",
    test_perturbations={"my_MEKi": ["trametinib", "selumetinib"]},
    panels="all_validated",   # tests against all 9 anchor (class, module, sign) triples
    n_perm=10000,
    seed=42,
)
print(result[["panel_class", "panel_module", "observed_z", "p_one_sided", "q_BH", "verdict"]])
```

Background null = all drugs in `adata.obs[pert_col]` not in `test_perturbations`. Need ≥30 background drugs (≈1000 in the validated Paper 4.1 procedure).

---

## Where to go next

- **`docs/EIGENSPACE_EVOLUTION.md`** — relationship between v0.1 12-PC and v0.3 9-PC eigenspaces. Read if comparing tool output to Paper 3.
- **`docs/SCOPE.md`** — explicit boundaries of the framework. What it answers, what it doesn't.
- **`README.md`** — full API reference and bundled data inventory.
- **`CHANGELOG.md`** — version history.
- **GitHub issues** — feature requests and bug reports.

## Common pitfalls

1. **Forgot to log-normalize before scoring.** `score_modules` expects `adata.X` to be log-normalized (`sc.pp.normalize_total(target_sum=1e4); sc.pp.log1p(adata)`). Raw counts will give wrong scores.

2. **Z-scoring across datasets.** If you score dataset A with `mean_zscore` and dataset B with `mean_zscore` separately and compare deltas, you get zero by construction. Use `mean_raw` for cross-dataset.

3. **Skipping validate_perturbation.** Module-level effects from a perturbation analysis are unreliable until you've checked random_200 + housekeeping + cell_cycle nulls. The Paper 4.5 v1.2 amendment exists exactly because skipping this check produced artifact-driven results.

4. **Reading capacity_floor as absolute prediction.** `saturated_blocked_up` means the cell cannot ramp the module UP under the operation's stimulus. It does NOT mean the module level can't change at all — specific signaling perturbations (e.g. atRA → UPR-ATF6) can drive saturated modules DOWN. The predictor is upward-asymmetric.

5. **Using `compare_to_references` for drug-disease matching.** It doesn't include drugs in v0.2 (cosine matching falsified by Paper 4.1). Use `pct.activity_layer_screen` for the validated drug operation.

6. **Trusting v0.1 PC2-PC3 interpretation in v0.3.** PC1 is preserved across versions. PC2-PC3 underwent rotation/reordering — see `EIGENSPACE_EVOLUTION.md`.
