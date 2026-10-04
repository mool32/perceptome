# Gate 0 — two layers on independent curations

## Setup

- instruction v1 (`AGENT_TASK_B_two_layers_gate0.md`) sha256 `9272659ea935c4429f1c15658c6d9b752f9e4348fb4fc5bb24abfc20dad061c5`;
  **instruction v2** (`AGENT_TASK_AB_continue.md`, §B replaces §3–4) sha256 `236bad96aa146192ba2101554d68f2025b3aa070134e51ee8a2b62b54dcd2ad6`
- Paper 7 project: `…/research/perceptual_modules/paper7_evolution` (not a git repository; no commit id). Pipeline scripts `step01…step05`, `step09` used as found.
- versions: python 3.12.0, numpy 1.26.4, scipy 1.18.1 (system interpreter, the one the pipeline runs under)
- parameters unchanged: 15 species, K = 0.6, P5 null = per-species coverage shuffle, 10⁴ permutations, SEED 42 / 43 as in the original code
- **Compara cache.** Each sandbox was seeded with the existing `symbol_resolution_cache.json` (203 entries) and `ortholog_cache.json` (214 entries). New symbols requested from Ensembl REST per catalog: C0 **0** (step02/03 finished in 0 s — the cache covered all 203 genes), C2 **158** (cache 203 → 361 symbols, ortholog cache 214 → 367; step02 428 s, step03 794 s), C3 **148** (→ 351 / 356; 322 s / 588 s), C4 **134** (→ 337 / 341; 456 s / 538 s). Ensembl REST was available throughout; no step failed.
- Catalog files found: `cancer_paper_release/data_specs/module_catalog.json` (v0.2, 43 modules, 203 core genes) and `perceptome_validation/registry/T5_catalog_{C2,C3,C4}.json`.
  **There is no separate `C1` file**; the catalog that matches the Paper 7 original is the project's own `module_catalog.json` — it is identical gene-for-gene to the current v0.3 catalog on all 43 shared modules (v0.3 only adds NPAS4). It is run here as **C0**.

### How the pipeline was pointed at another catalog

Every Paper 7 script derives all its paths from a single constant `BASE`, and reads the catalog from
`BASE/cancer_paper_release/data_specs/module_catalog.json`. So each catalog got a sandbox mirroring that layout, and the
copies `*_gate0.py` differ from the originals in **exactly one line each — the `BASE` assignment** (`gate0_scripts.diff`,
24 lines total for all six scripts). Nothing else was touched: no parameter, no statistic, no threshold. Output goes to the
sandbox, so the original Paper 7 data files were not overwritten.

Note that the hard-coded `BASE` in the originals is `/Users/teo/Desktop/research/perceptual_modules`, a path that **no longer
exists** (the project has moved), so the scripts cannot run at all without this change — it is both the sanctioned catalog-path
change and a precondition for running.

## Reproduction (C0)

| module | z original | z rerun | \|Δz\| | age original | age rerun |
|---|---|---|---|---|---|
| HSF1 | 3.875 | 3.886 | 0.011 | 1500 | 1500 |
| ERK/MAPK | 2.877 | 2.863 | 0.014 | 1500 | 1500 |
| PI3K/PTEN | 2.869 | 2.885 | 0.017 | 1500 | 1500 |
| AMPK | 2.878 | 2.897 | 0.019 | 1500 | 1500 |
| p53 | −2.133 | −2.146 | 0.013 | 684 | 684 |
| **VDR** | −2.163 | −2.174 | 0.011 | **435** | **797** |
| cGAS-STING | −3.109 | −3.146 | 0.038 | 435 | 435 |
| NFAT | −3.126 | −3.168 | 0.042 | 435 | 435 |

Over all 43 modules: **max \|Δz\| = 0.050**, median 0.014 — well inside the ±0.3 tolerance.

**The age mismatch and what causes it.** Core-gene lists, `symbol_to_ensembl.csv`, `ortholog_cache.json` (214 entries) and the
gene-level ortholog matrix are byte-identical between the stored run and the rerun, yet the stored `module_coverage_matrix.csv`
has VDR exactly one gene lower in *every* species (mouse 0.8 vs 1.0, chicken 0.4 vs 0.6, fly 0.4 vs 0.6). Recomputing VDR
coverage from the current cache gives 5/5 in mouse and 3/5 in chicken and fly, i.e. the rerun's values. So the stored
`module_coverage_matrix.csv` / `module_presence_matrix.csv` / `module_ages.csv` (and the ages inside `robust_cohorts.json`) were
produced with an **earlier state of the ortholog cache**, which has since been extended. Four of 43 modules cross the K = 0.6
threshold differently as a result, all of them in the direction of greater age:

| module | stored | rerun |
|---|---|---|
| NRF2 | 435 (D.rerio, 6 species) | 1500 (A.thaliana, 8) |
| PXR/CAR | 797 (C.elegans, 6) | 1500 (D.discoideum, 10) |
| Type I IFN | 684 (C.intestinalis, 6) | 1500 (D.discoideum, 8) |
| VDR | 435 (D.rerio, 5) | 797 (D.melanogaster, 7) |

Of the Tier-1 eight, only VDR is affected. The run was stopped at this point and the discrepancy reported; **the user chose
option 1**: reproduction accepted on z, ages taken as recomputed on the current cache, continue to C2–C4. All numbers below
therefore use the current cache for every catalog including C0, so the four catalogs are treated identically.

## Per catalog

Full per-module tables: `per_module_C0.tsv`, `per_module_C2.tsv`, `per_module_C3.tsv`, `per_module_C4.tsv`
(columns: module, n_core, age_mya, z, observed_intact, null_mean, null_std, category).

| catalog | modules | G0-1: n ancient (≥1500) | n young (≤800) | median z ancient | median z young | one-sided MW p | G0-1 |
|---|---|---|---|---|---|---|---|
| C0 (Paper 7) | 43 | 18 | 23 | +1.87 | −1.14 | 1.4 × 10⁻⁶ | pass |
| **C2** | 49 | 12 | 35 | +2.57 | −0.81 | 2.7 × 10⁻⁷ | **pass** |
| **C3** | 48 | 12 | 32 | +3.32 | −0.40 | 8.6 × 10⁻⁷ | **pass** |
| **C4** | 47 | 10 | 36 | +2.13 | −0.00 | 1.0 × 10⁻⁶ | **pass** |

Every catalog has ≥ 5 modules in both groups, so none is n/a.

### G0-2 (descriptive only, does not affect the verdict)

Candidates are modules of Ck whose T5 `sensor`/`tf` fields contain a gene of the Paper 7 module's `sensor_genes` ∪ `tf_genes`
(`modules_v03.json`; `core_genes` where no `sensor_genes` exists); among candidates, the highest core-gene Jaccard.

| Tier 1 | C2 → module (J, z, age) | C3 | C4 |
|---|---|---|---|
| HSF1 | Heat shock/proteotoxic (HSP90AA1-HSF1) — 0.571, **z +5.14**, 1500 | Heat shock/cytosolic misfolding (HSF1) — 0.429, **+3.96**, 1500 | Cytosolic proteotoxic stress / heat shock (HSF1) — 0.667, **+4.37**, 1500 |
| ERK/MAPK | EGF/ErbB growth factor (EGFR-ELK1) — 0.100, −0.13, 797 | EGF/EGFR-ERK — 0.100, +0.65, 797 | EGF-family growth factor (ERK-AP-1) — 0.222, +0.89, 797 |
| PI3K/PTEN | not matched | Energy/AMP stress (AMPK-FOXO3) — 0.111, +3.90, 1500 | not matched |
| AMPK | not matched | not matched | not matched |
| p53 | DNA double-strand break (ATM-TP53) — 0.375, −0.15, 797 | DNA double-strand breaks (ATM-p53) — 0.250, −0.06, 1500 | DNA double-strand break (ATM-p53) — 0.222, +0.89, 1500 |
| VDR | Vitamin D (VDR) — 0.429, **−1.80**, 797 | Vitamin D (VDR) — 0.500, **−1.05**, 797 | Vitamin D (VDR-CYP24A1) — 0.500, **−0.81**, 797 |
| cGAS-STING | Cytosolic DNA STING (CGAS-IRF3) — 0.111, −0.81, 797 | Cytosolic DNA (cGAS-STING-IRF3) — 0.111, −1.07, 797 | Cytosolic dsDNA sensing (cGAS-STING-IRF3) — 0.111, −0.88, 797 |
| NFAT | Gq-calcium NFAT (AGTR1-NFATC1) — 0.375, +1.95, 1500 | Membrane stretch (PIEZO1-calcineurin-NFAT) — 0.375, +1.73, 1500 | Gq-coupled GPCR / calcium (calcineurin-NFAT) — 0.375, +1.98, 1500 |

Signs of z: HSF1 is positive in all three (and the largest z in each catalog); VDR and cGAS-STING are negative in all three;
ERK/MAPK, p53 and NFAT change sign between catalogs; AMPK is never matched by this rule, PI3K/PTEN only in C3 (and to an
AMPK-named module, Jaccard 0.111 — a weak match). The mapping is as fragile as the earlier name check suggested; the table is
provided for manual inspection as instructed, and nothing in the verdict rests on it.

## Verdict

**PASS** — G0-1 passes in **3 of 3** alternative catalogs (C2, C3, C4); none is n/a.

Secondary (descriptive), Spearman ρ(age, z) over all modules of the catalog:

| catalog | ρ | p |
|---|---|---|
| C0 | +0.797 | 1.5 × 10⁻¹⁰ |
| C2 | +0.851 | 1.0 × 10⁻¹⁴ |
| C3 | +0.862 | 3.3 × 10⁻¹⁵ |
| C4 | +0.876 | 8.2 × 10⁻¹⁶ |

The `B_interoceptive` comparison is dropped per instruction v2 (T5 modules carry no such category).

## Deviations

1. **`BASE` changed in the script copies** — one line per script, six scripts per catalog (`gate0_scripts.diff`). Required twice
   over: it is the sanctioned catalog-path change, and the original hard-coded path no longer exists.
2. **T5 catalogs adapted to the input format** step01 expects: each module's `genes` written as `core_genes`, and
   `category` set to `"unknown"` because the T5 curations have no category field. Nothing in the statistics uses `category` —
   step01 only copies it into `module_genes.csv`, step05 only reports it.
3. **C0 reproduction accepted on z, with recomputed ages** — the user's decision (option 1) after the stop described above.
   Consequence to keep in view: C0's own age groups here are not the ones in `robust_cohorts.json` (4 of 43 modules moved).
4. Ages, z and p for C0 in this report come from the rerun, not from the stored files, so that all four catalogs are computed
   under identical conditions.

Nothing was committed or pushed; all outputs are under `paper7_evolution/results/gate0_two_layers/`.

## Confirmation

- pipeline unmodified except the catalog path (diff attached): **yes** — `gate0_scripts.diff`, 24 lines for all six scripts, every
  hunk the `BASE` assignment and nothing else
- no additional analyses: **yes** — only the criteria of instruction v2 (G0-1, the descriptive G0-2 mapping, Spearman ρ) and the
  reproduction comparison the instruction asks for. The VDR diagnosis was tracing an input mismatch the instruction required me to
  report, not a new analysis. No plots.

## Attached

`per_module_{C0,C2,C3,C4}.tsv` · `gate0_summary.json` · `gate0_scripts.diff` · `run_{C0,C2,C3,C4}.log`, `run_C2C3C4.out` ·
`setup_gate0.py`, `evaluate_gate0.py` · per-catalog sandboxes `C*/paper7_evolution/{data,results}` with the full pipeline outputs.
