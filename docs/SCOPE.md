# Scope — what perceptome answers, what it doesn't

This page draws the boundary explicitly. The framework has been pre-registered, tested, and (in many cases) falsified across multiple substrate domains. Knowing the boundary is part of using the tool well.

The discipline: *every framework claim with broad scope was tested. Falsified and retracted claims are recorded with the same prominence as validated ones.* Five framings of the drug-pharmacology question were tested in a single day and all killed; one substrate-series predicted scope condition was empirically narrowed; one cancer-attractor formulation split into a static-killed / dynamic-passed verdict; and the convergence paper (Paper 3) retracts its own earlier baseline-gating reading and withdraws several arguments (section 1).

This page exists so that **users do not waste effort on operations the framework does not support**.

---

## 1. Convergence layer (Paper 3)

Source: Spiro (2026), *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system* ([Zenodo 10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130)). Implemented by `pct.convergence`; the numbers are in `pct.load_convergence_reference()`; the user guide is [`CONVERGENCE.md`](CONVERGENCE.md).

### Validated

- **Malignant cells of 25 cancers occupy a restricted region of the fixed normal-cell space.** Mean pairwise cosine distance 0.16 (0.166) among malignant states against 0.99 among their near-orthogonal normal origins. The region survives leave-one-cohort-out (core fraction 25.4–28.4% vs 27.1%; the cytotrophoblast's recurrence rank is 2 in all 25 panels), dropping five of the 44 modules (17.5 ± 2.5%, never below 11.7%), dimensionality (six vs nine dimensions), the assignment rule (soft, top-k) and a 36-way specification curve under the single-cell calibration.
- **Most of that contraction is shared transcriptional output; a smaller module-shape convergence remains.** The raw contraction is matched by a matched-averaging null (0.166 vs 0.129). With magnitude removed the malignant centroids cluster more tightly than the null in two of three scalings (raw module space: 0.33 vs 0.49, below all 200 draws; reference scaling: 0.026 vs 0.034, P = 0.03); under a third scaling the advantage disappears (0.40 vs 0.38). The paper reports all three.
- **The placental cytotrophoblast is the lineage-independent anchor under this catalog.** A top target in 18 of 25 cancers; 11.1% of malignant cells vs 0.66% by chance (16.9-fold); the cytotrophoblast's assignment rises with transcriptome complexity (6.8% → 21.6%) while the megakaryocyte's falls (13.2% → 3.5%), tracking the quarantined erythrocyte — so the megakaryocyte is carried as a secondary target only. The region is not the manifold's centre (the medoid receives zero convergence) and not the germline.
- **Within patients, malignant epithelium occupies the core more than matched normal epithelium** — the paper's most defensible claim. Colorectal 18.9% vs 3.2%, 43 of 46 patients (+10.5 points with magnitude removed, 44/46; +12.9 points resampled to matched complexity, 43/46); an independent colorectal dataset 23/25; kidney 6/7 and liver 5/8 same-signed; pooled over four datasets 77 of 86 patients (+13.0 points, mixed model with a dataset random intercept); held-out gastric 18/22; lung 10/10. Not specific to this catalog: KEGG, Reactome, Hallmark and PROGENy reproduce it at the same size or larger.
- **The state is physiological, not cancer-specific.** Healthy tissue occupies it to a tissue-dependent degree (lung epithelium 16.1%, neural 7.6%, kidney 0.5%), and fibrotic lung — the strictest non-malignant comparator — does not enter it (+3.8 points raw, +1.3 shape-only; n.s.), whereas malignant lung does (+15.7 points).
- **The target-defining signature is reproduced by cancer.** The programs that make the core unusual among normal cells (about a dozen suppressed brakes: Circadian, Type I IFN, JAK-STAT, PI3K/PTEN, Notch, Wnt, SREBP …) correlate with what cancers contribute as they converge: module-level Spearman ρ = +0.88 (bootstrap 95% CI 0.78–0.93), gene-level ρ = 0.56. Type I IFN is the one brake cancer does not reproduce.
- **Malignant cells already in the state add stress and growth programs** rather than moving further along the defining axis (ρ = +0.36; colorectal only).

### Retracted or superseded in the paper

- **Baseline gating** — the reading that the within-patient increase is present where normal tissue sits low and absent in high-baseline lung — and with it the prostate "saturated high-baseline" justification. The apparent absence in lung was a cell-selection artifact (the Census lung "tumour" arm was ~9% annotated-malignant cells and ~91% normal-adjacent epithelium). No baseline-gating claim is made.
- **The fold-change framing** ("about sixfold"): inflated by a small normal denominator. The effect is reported in percentage points (~ +10, flat against baseline on the additive scale).
- **The comparison against random 25-cell-type draws** (a Stayton-type test) and the fold-contraction it produced; the null ignored the averaging itself.
- **The shared-polyploidy argument** for the two core cells (the anchor is the mononuclear diploid villous cytotrophoblast).
- **"Malignant occupants run cGAS-STING hotter"**: the module was scored on a three-gene fragment lacking both sensors (Note S-Q).
- **The argument that malignant-cell contamination makes the reported numbers a lower bound** (tested on HCC and BRCA; not borne out).
- Four candidate organizing principles of the state — driver-gene enrichment, module evolutionary age, an exteroception–interoception axis, immune privilege — were tested and falsified.

### Suggestive or open

- **Time resolution.** Lung-adenocarcinoma cells under erlotinib move into the state and do not return (shape-only core 29% → 53% → 30% → ~40%); a cell line without lineage tracing, so movement cannot be separated from selection. The coordinate is not moved by single-gene knockdown (Perturb-seq, K562, an unfavourable substrate; combinatorial knockdowns untested).
- **Pre-invasive prediction.** Holds in colorectum (adenoma above matched normal in 3/3) and directionally in breast (DCIS between normal and invasive); fails in esophagus (non-dysplastic Barrett's below normal squamous in 4/4); dysplasia is untested.
- **Metastasis.** The matched lymph-node metastasis sits lower than the primary (breast 21.8% → 13.4%, 6/6) — nodal not distant metastases, small numbers.
- **Prostate** is an open negative: two paired cohorts disagree about the baseline itself, neither shows an increase.
- **Reprogramming.** The fraction of cells assigned to the cytotrophoblast rises tenfold late in mouse fibroblast-to-iPSC reprogramming (ρ = +0.76); the starting fibroblasts score anomalously high and the dataset has no cell-type labels, so the paper calls it suggestive only.

### Not supported

- **A proven attractor, or convergent evolution.** The test is a static geometric one: a convergent endpoint across origins, not a measured trajectory. The paper says "convergence, not convergent evolution".
- **That the cytotrophoblast is *the* convergence target independent of the catalog.** Three independent reconstructions of the catalog (three instances of one model family, sharing a prior) reproduce the region but name the gastric progenitor; the anchor is a statement about this catalog, and the 154-type gallery is sparse in progenitor and fetal cells.
- **Inheritance of immune privilege.** Across 62 colorectal tumours the core fraction is unrelated to T/NK infiltration (ρ = +0.04 shape-only, +0.13 raw).
- **Ranking cancers by convergence fraction.** The ordering is sensitive to malignant-label quality (one colon dataset moved from 13% to 36% after an author-validated relabeling) and to platform effects.
- **A cancer-specific destination.** Normal tissues occupy the state; the unit of the result is a physiological state whose occupancy malignancy increases.

### Preconditions for using the assay

- **Per-cell data need the single-cell-scaled calibration** (bundled; the default of the per-cell pipeline). Under the pseudobulk calibration cells collapse onto low-complexity targets (core fraction 7% vs 17%; 99.6% of independent normal single cells onto erythroid targets).
- **≥ 1,000 detected genes per cell**, the erythrocyte and platelet quarantined, **≥ 20 cells per arm per patient**, the patient as the unit of replication.
- **A matched, correct comparator** (e.g. crypt-base progenitors for colorectal cancer) and **malignant-cell purity** (author-annotated compartments; CNV-inference floors across all 25 cohorts are deferred in the paper).
- **Gene symbols**: three catalog genes have been renamed; the corrected catalog has not been re-derived end to end (Note S-Q).
- **The criteria are curation records.** The four operational criteria were not validated empirically as a selection rule (feedback genes are not enriched for cancer-driver status); a feedback element is named for 38 of 44 modules and, for ten of them, is not among the scored genes. The curated gene list was assembled with cancer biology in view: reshuffling the same 206 genes into random modules still gives 14% convergence against 6% for random genome-wide panels and 21% for the real modules, so roughly two-thirds of the effect is carried by the gene set, not by its grouping.

---

## 2. Geometry layer (eigenspace)

### Validated

- **Project any cell into the fixed 9-PC reference space built from 154 normal HPA cell types** (standardized coordinates: module scores z-scored against the reference calibration, then the loadings). Nine components by the Kaiser criterion (73.8% of variance); PC1–PC6 bootstrap-stable (0.63–0.99), PC7–PC9 probable (0.40–0.54). The axes are individually interpretable: PC1 runs from inflammatory cells to terminal germline cells, PC2 from quiescent to proliferative cells, PC3 from signaling-dominated immune cells to metabolic epithelia. The first Paper 3 draft's labels — PC1 as "perception breadth", PC4 as a "cancer convergence axis" — are not used by the final paper, which describes convergence as a high-dimensional nearest-neighbour property rather than an axis or a 2-D location.
- **Calibration matters.** Both calibrations are bundled: `pseudobulk` for cell-type-level scores, `single_cell_scaled` for per-cell data (see section 1).
- **Compute (tumor − normal) shift in eigenspace coordinates** (`pct.project_shift`) and compare to known directions (legacy cancer capacity-direction — Paper 4.2 P3 PASS in Sun 2021 paired HCC; the 8-cell attractor cluster of that analysis is superseded by section 1).
- **Per-cell or per-cluster scoring**, with `mean_raw` (default), `mean_zscore` (within-dataset only), or `scanpy_corrected` (proliferation contexts).

### Recomputed for v0.3

- Disease (RA/AD/IPF/DKD) and aging (inflammaging/collapse) reference vectors rebuilt from raw scRNA-seq on v0.3 eigenspace and expressed in the standardized coordinates. See `CHANGELOG.md` for details + biological refinement (shared aging modules narrowed from 4 to 2). The paper excludes its rheumatoid-arthritis atlas from the disease controls because its median is 286 detected genes per cell, below the 1,000-gene floor: check the complexity of the RA data behind the RA vector before relying on it.
- Cross-version comparison: PC1 ports across v0.1↔v0.3 with high fidelity; PC2-PC3 underwent rotation (see `EIGENSPACE_EVOLUTION.md`).

### Not supported

- **Drug-disease cosine matching in eigenspace.** Paper 4.1 falsified in three independent formulations:
  - v1.0 cosine of 12D vectors: 23/50 pairs concordant (chance), p=0.76
  - v1.1 stress vs specificity subspace decomposition: no subspace recovers signal
  - v1.2 module-level overlap test: 27/50 pairs concordant, hypergeometric q<0.10 in 0/50
  - The `compare_to_references()` API does not include drugs. Use `pct.activity_layer_screen` for the validated drug operation.

---

## 3. Perceptivity layer (capacity)

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

- 44 modules cover general signaling. Tissue-specific modules (HNF4α liver, MYOD muscle) not yet in catalog (planned). NPAS4 added as 44th canonical because it passes all four operational criteria; HPA shows it neuron-enriched (mean A_NPAS4 in neurons = 2.97 vs ≤2.13 in other classes). It does not shape the cross-tissue space: a 43-module rebuild reproduces all nine components with eigenvector cosine ≥ 0.95 (paper Note S2).
- The four criteria record how the modules were curated; the paper could not validate the feedback criterion empirically, and the curated gene list carries a residual cancer-biology prior (section 1).

---

## 4. Validity layer

### Validated

- **Random-200 / housekeeping / cell-cycle scorecard** for any perturbation-vs-control comparison. Standard practice since Paper 4.5 v1.2 amendment caught a baseline-shift artifact in the adult-intestinal-epithelium analysis (random-200 d=−0.60 → ARTIFACT FAIL → switched to scanpy_corrected scoring).
- **Switching to `scanpy_corrected`** when random_200 fails. This is the documented fix for the proliferation-baseline confound.

### Not supported

- The validity scorecard does not catch all confounds. Cross-platform batch effects, cell composition shifts within clusters, and ambient RNA contamination are not detected. Use upstream tools (scanpy QC, cellxgene-census batch correction) for those.

---

## 5. Drug perturbation analysis

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

## 6. Methodological commitments

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
