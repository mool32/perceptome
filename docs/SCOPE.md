# Scope — what perceptome answers, what it doesn't

This page draws the boundary explicitly. The framework has been pre-registered, tested, and (in many cases) falsified across multiple substrate domains. Knowing the boundary is part of using the tool well.

The discipline: *every framework claim with broad scope was tested. Falsified claims are recorded with the same prominence as validated ones.* Five framings of the drug-pharmacology question were tested in a single day and all killed; one substrate-series predicted scope condition was empirically narrowed; one cancer-attractor formulation split into a static-killed / dynamic-passed verdict.

This page exists so that **users do not waste effort on operations the framework does not support**.

---

## 1. Geometry layer (eigenspace)

### Validated

- **Project any cell into 9-PC HPA-derived perceptual space.** PC1 = perception breadth (Paper 3 PC1 preserved, |cos|=0.90). PC4 carries the cancer convergence direction (Paper 3 PC4 preserved, |cos|=0.71). PC1-PC6 bootstrap-stable.
- **Compute (tumor − normal) shift in eigenspace coordinates** and compare to known directions (cancer attractor — Paper 4.2 P3 PASS in Sun 2021 paired HCC).
- **Per-cell or per-cluster scoring**, with `mean_raw` (default), `mean_zscore` (within-dataset only), or `scanpy_corrected` (proliferation contexts).

### Recomputed for v0.3

- Disease (RA/AD/IPF/DKD) and aging (inflammaging/collapse) reference vectors rebuilt from raw scRNA-seq on v0.3 eigenspace. See `CHANGELOG.md` for details + biological refinement (shared aging modules narrowed from 4 to 2).
- Cross-version comparison: PC1 ports across v0.1↔v0.3 with high fidelity; PC2-PC3 underwent rotation (see `EIGENSPACE_EVOLUTION.md`).

### Not supported

- **Drug-disease cosine matching in eigenspace.** Paper 4.1 falsified in three independent formulations:
  - v1.0 cosine of 12D vectors: 23/50 pairs concordant (chance), p=0.76
  - v1.1 stress vs specificity subspace decomposition: no subspace recovers signal
  - v1.2 module-level overlap test: 27/50 pairs concordant, hypergeometric q<0.10 in 0/50
  - The `compare_to_references()` API does not include drugs. Use `pct.activity_layer_screen` for the validated drug operation.

---

## 2. Perceptivity layer (capacity)

### Validated

- **Capacity floor predictor (upward-asymmetric).** A_baseline > 4.5 reliably predicts a cell cannot ramp the module UP further. Closed across **three independent paradigms** in 4 datasets:
  - Paper 4.5 Elmentaite Gut Cell Atlas (in vivo adult homeostatic intestinal epithelium): Goblet UPR-ATF6 d=−0.10 ✓
  - Paper 4.6 Mead 2022 GSE148524 (cocktail differentiation): max log2FC=0.234 ✓
  - Paper 4.7 Lukonin 2020 GSE147135 (regeneration paradigm): UPR-ATF6 |log2FC|=0.023, NRF2 |log2FC|=0.030 ✓
  - HPA static atlas (reference, 154 cell types) ✓
- **Two-factor framework Factor 1.** The tool computes `capacity` from HPA. The full prediction is `engagement = capacity × operation_intensity`; Factor 2 (operation intensity) is context-specific and not in the tool.

### Refined scope

- **The capacity-floor predictor is upward-asymmetric, not absolute.** Paper 4.8 holdout PARTIAL (3/5 PRIMARY) revealed that specific signaling perturbations (atRA via PPAR cross-talk) can DOWNregulate saturated modules (UPR-ATF6 −0.335, UPR-PERK −0.247). This refines the predictor: cells with saturated baselines cannot ramp UP, but they can be actively suppressed DOWN by specific perturbations. Recorded as scope refinement, not refutation.

### Empirically narrowed

- **Architecture engagement requires an active context.** Paper 4.5 epithelial closed FAIL (1/8 PRIMARY): the operation/perception/infrastructure architecture is observed in:
  - Active differentiation (neurons LTP, plasmablast tonsil, Th1)
  - Acute remodeling (muscle hypertrophy)
  - Regeneration (organoid)
  - NOT in steady-state homeostatic terminally-differentiated cells (adult intestinal goblet, Paneth, enterocyte)
- Tool provides the capacity prediction; whether your cells are in an active context is a biological question the user must answer.

### Catalog scope limitation

- 44 modules cover general signaling. Tissue-specific modules (HNF4α liver, MYOD muscle) not yet in catalog (planned for v0.3). NPAS4 added as 44th canonical because it passes all four perceptual criteria; HPA shows it neuron-enriched (mean A_NPAS4 in neurons = 2.97 vs ≤2.13 in other classes).

---

## 3. Validity layer

### Validated

- **Random-200 / housekeeping / cell-cycle scorecard** for any perturbation-vs-control comparison. Standard practice since Paper 4.5 v1.2 amendment caught a baseline-shift artifact in the adult-intestinal-epithelium analysis (random-200 d=−0.60 → ARTIFACT FAIL → switched to scanpy_corrected scoring).
- **Switching to `scanpy_corrected`** when random_200 fails. This is the documented fix for the proliferation-baseline confound.

### Not supported

- The validity scorecard does not catch all confounds. Cross-platform batch effects, cell composition shifts within clusters, and ambient RNA contamination are not detected. Use upstream tools (scanpy QC, cellxgene-census batch correction) for those.

---

## 4. Drug perturbation analysis

Paper 4.1 closed 2026-05-09. The full audit trail is in `PAPER4_1_CANONICAL_RESULTS.md` and `killed_hypotheses.csv` in the source tree.

### Validated (6 + 3 anchors)

- **6 mechanism-pathway recoveries** via activity-layer scoring against background null:
  - MEKi → ERK/MAPK ↓ (Block 5 z=−1.026, q=0.0002; Block 2 holdout PASS, drugs PD-98059, mek1-2-inhibitor)
  - Proteasomei → UPR-PERK ↑ (Block 5 z=+1.209, q=0.0002; Block 2 holdout PASS, ixazomib)
  - CDKi → Cell Cycle ↓ (Block 4 LD=+3.30; not in Block 5)
  - EGFRi → ERK/MAPK ↓ (Block 4 LD=+1.07)
  - PI3Ki → PI3K/PTEN ↑ FOXO targets (Block 4 LD=+0.81)
  - IKKi → NF-κB ↓ (Block 4 LD=+0.79)
- **3 stress-axis positive controls** (Block 5 v1.2; assayed for null calibration, not claimed as primary findings):
  - HSP90i → HSF1 ↑ (z=+2.149, q=0.0002)
  - HSP90i → UPR-ATF6 ↑ (z=+0.736, q=0.0004)
  - HSP90i → NRF2 ↑ (z=+0.510, q=0.038)
- These are the rows of `pct.drug_anchors()`.

### Falsified (11 framings)

The 11 pre-registered falsifications cover the obvious-looking drug operations that **don't work**, even though many of them sound like exactly what users would want:

1. **Drug-disease cosine matching in eigenspace** (3 formulations: v1.0 cosine, v1.1 subspace decomposition, v1.2 module overlap)
2. **Drug-class mechanism recovery** via panel geometry at 7/10 threshold (5/10 PASS)
3. **Within-panel-shuffle null** (statistically miscalibrated — replaced by background null v1.2)
4. **TD/PT layer-differential** as universal class predictor (6/18 correctly placed)
5. **Readiness-layer timescale rescue** (24h vs 6h paired Δ): wrong direction on HDACi
6. **TF-autoregulation as 1st-class layer**: K=1/8 classes, binomial P(K≥1)=0.57 = noise
7. **Equilibrium-instrument linear-response framing**: K=1/6 classes, Proteasomei wrong direction
8. **"Clean perceptome signature predicts FDA approval beyond selectivity"**: ΔR²=+0.002, interaction wrong direction
9. **"Snapshot perceptome blind to pulsatile dynamics"**: 1/6 rescues actually pulsatile, conjunction criterion mathematically infeasible

### Not supported

- Any of the 11 falsified operations above. Implementing them in the tool would invite users to do exactly the analysis the data already refused.
- Mechanism deconvolution from drug centroid geometry. Uses panel-design consequence, not deeper mechanism.

### Five interpretive framings tested in one day, all killed

Pattern noted explicitly in Paper 4.1 §4: every theoretical interpretation of WHY the 6 surviving rescues survive (TD/PT mechanism, timescale, TF dynamics, thermodynamic regime, program engagement, snapshot dynamics-blindness) has been tested and falsified at noise floor. **The honest framing is panel-design consequence, no deeper unifying theory.** The tool's drug subpackage implements only the descriptive operation, not any of the killed interpretations.

---

## 5. Methodological commitments

The tool reflects these disciplines from the underlying research:

1. **Pre-registration before data inspection.** Threshold and verdict criteria locked in pre-reg files (SHA-archived) before computation. See `paper4.x/preregistration/` directories in the source tree.
2. **Asymmetric thresholds.** Positive results require multiple confirmation steps (random-null + holdout + scramble). Negative results are self-correcting and reported once.
3. **Scramble null where applicable.** `activity_layer_screen` uses background null + scramble-equivalent permutation.
4. **Failed hypotheses recorded.** `killed_hypotheses.csv` in the research tree documents 21 falsifications across the program. The drug subpackage in the tool implements the 9 surviving anchors and explicitly does not implement the falsified operations.
5. **No post-hoc rescue.** When a pre-reg fails, it is recorded killed with full numbers, not reframed.

These disciplines are why the tool is narrower than it could be. Releases without this discipline tend to grow features that do not survive replication.

---

## When to NOT use perceptome

- **Single-gene analyses.** Use scanpy.tl.score_genes or DESeq2 directly.
- **Cell type identification or label transfer.** Use scanpy/CellTypist/Symphony.
- **Trajectory inference.** Use Slingshot / scFates / palantir.
- **Drug repositioning via geometric matching.** Falsified by Paper 4.1 — see drug section above.
- **Phosphorylation cascade readout.** Perceptome is a transcriptional instrument; mTOR/MEK direct phosphorylation outputs are SC3 (scope condition 3 in Paper 4.1) — invisible to transcriptome.
- **Clinical decision support.** This is a research tool. No clinical validation has been performed.
