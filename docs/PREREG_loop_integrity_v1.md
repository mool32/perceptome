# Pre-registration v1: does a broken feedback loop leave a noise signature?

**Status:** LOCKED (see `docs/PREREG_LOCK.txt`, git tag `prereg-loop-v1`). Locked before any
expression value of the screen was read; only `obs`/`var` metadata had been inspected.

**The analysis is the code.** `scripts/14_loop_integrity_prereg_v1.py` implements every rule
below and is locked with this file. Where prose and code disagree, the code governs and the
disagreement is reported as a deviation.

## 1. Question

Negative feedback compresses cell-to-cell variability of a pathway's output (Becskei & Serrano
2000). If that holds and is readable in single-cell RNA, then removing a module's feedback
("break") should raise the dispersion of its output genes **beyond what the mean change
explains**, while removing its input ("input_down") should not. Pathway-activity methods
see only the mean; this test asks whether a second, mechanistic readout exists.

## 2. Data

- Replogle et al. 2022, K562 genome-wide CRISPRi Perturb-seq (`K562_gwps_raw_singlecell_01.h5ad`,
  figshare article 20029387 v1). Raw UMI counts.
- Access: full download of the file; verify the figshare md5 before use and record sha256.
- Metadata already extracted (pseudobulk obs/var only): `data/evidence/`.

## 3. Analysis set (frozen from metadata only)

Eligibility of a knockdown (KD): >= 50 filtered cells **and** author `fold_expr` <= 0.50
(>= 50 % knockdown). Eligibility of a loop: >= 1 eligible KD in `break` (or `input_up`) and in
`input_down`, and >= 3 output genes (`pct.loop_output_genes`) measured in the K562 readout.

| loop | break KDs | input_down KDs | input_up KDs |
|---|---|---|---|
| ERK/MAPK | SPRY2, SPRY4, SPRED1 | SOS1, RAF1, MAP2K1, MAPK1 | NF1 |
| HIF | VHL, EGLN1 | HIF1A, ARNT | — |
| JAK-STAT | SOCS1, CISH, PTPN2 | STAT5A, STAT5B, JAK2, JAK1 | — |
| SREBP | INSIG1, INSIG2 | SCAP, SREBF2, MBTPS1, MBTPS2 | — |
| UPR-IRE1 | DNAJB9 | ERN1, XBP1 | — |
| UPR-PERK | PPP1R15A | EIF2AK3, ATF4 | — |

Primary: 12 break KDs, 18 input_down KDs (6 loops).
Excluded by the eligibility rule (reported, not analysed in the primary): DUSP6 (absent),
SOCS3 (40 cells), GRB2 (fold 0.77), STAT3 (0.63), SREBF1 (0.89), HSPA5 (0.57, so UPR-ATF6 has
no eligible break), HSF1 (0.84, so HSF1 has no eligible input_down), NFE2L2 (0.63, so NRF2 has
no eligible input_down), TSC2 (1.12).

Secondary sets:
- **break-only:** NRF2 (KEAP1, CUL3), HSF1 (HSP90AB1, HSP90AA1).
- **input_up vs input_down:** mTOR (TSC1 vs RPTOR, MTOR, RHEB); ERK (NF1).
- **negative control:** p53 (MDM2, PPM1D, MDM4, TP53, ATM). K562 carries a frameshift insertion in
  TP53 and expresses no functional p53 (Law et al. 1993, Leuk Res 17:1045; Cellosaurus CVCL_0004).
  The run records the Cellosaurus TP53 annotation; if it does not report a loss-of-function
  TP53 variant, N1 is reported but not used in the decision.
- **shared-node specificity:** MBTPS1/2 (SREBP and UPR-ATF6), ATF6 (UPR-ATF6 only), SCAP (SREBP only).

Output genes = `pct.loop_output_genes(module)` ∩ K562 readout, **excluding the knocked-down gene
itself**.

## 4. Cells and preprocessing

- Perturbed cells: all filtered cells of the KD (all guide groups pooled).
- Controls: non-targeting cells from the same gem groups as the perturbed cells, 5 per perturbed
  cell, sampled with seed 0; the same control cells are reused for every module of that KD.
- Depth: each cell's raw counts binomially thinned to a common depth D, where D is the 10th
  percentile of total UMI over all selected non-targeting cells (fixed once, before any
  perturbed cell is loaded); cells below D are dropped and counted.

## 5. Metrics (per KD k, module m)

- **M1 — mean shift.** For each output gene, log2((mean_k + 0.1) / (mean_c + 0.1)) on thinned
  counts; M1 = mean over output genes.
- **M2 — excess dispersion.** For each output gene with mean >= 0.1 in both groups,
  method-of-moments overdispersion phi = (var − mean) / mean²; d_g = ln((phi_k + 0.01) /
  (phi_c + 0.01)); M2 = median of d_g over output genes.
- **Z2 — specificity of M2.** 1000 random gene sets of the same size from readout genes outside
  the loop layer, matched to the output genes on control-mean decile; Z2 = (M2 − mean_null) /
  sd_null, same KD, same cells.
- **M3 — coupling (exploratory).** Spearman ρ across cells between the output score and the
  feedback score (mean z of log1p thinned counts; feedback set minus the KD gene); Δρ = ρ_k − ρ_c.
  Computed only where >= 1 measured feedback gene remains.
- Uncertainty: 500 bootstrap resamples of cells within each group, seed 0.

## 6. Gates (run first, in this order)

- **G1 — readout reaches the loop.** Per loop: median M1 over its input_down KDs < 0, and the
  bootstrap 95 % CI of M1 excludes 0 for at least one of them. Loops failing G1 leave the
  primary analysis and are reported as "loop not readable in K562".
- **G2 — metric is not an artifact.** 50 pseudo-perturbations, each 150 non-targeting cells
  versus other non-targeting cells from the same gem groups, scored on each primary loop's
  output genes: |Z2| > 1.96 in <= 10 % of them. **If G2 fails, stop:** the dispersion metric is
  dominated by sampling noise and nothing below is interpreted.
- **Kill switch:** fewer than 3 loops pass G1 → stop, report as "readout too weak at this depth".

## 7. Hypotheses

- **H1 (primary).** Over break KDs of loops passing G1: mean Z2 > 0 (sign-flip permutation,
  10,000, one-sided p < 0.05) **and** Z2(break) > Z2(input_down) (Mann–Whitney, one-sided
  p < 0.05).
- **H1-spec (required with H1).** Linear model over all primary KDs: Z2 ~ class + |M1|; the break
  coefficient > 0 with one-sided p < 0.05. Guards against "break raises the mean, and a higher
  mean alone shifts dispersion".
- **N1 (negative control, prediction).** In TP53-null K562, removing p53's feedback produces no
  noise signature: Z2 < 2.128 on p53 outputs for each break KD (MDM2, PPM1D, MDM4; one-sided,
  Bonferroni over 3). M1 and the input KDs (TP53, ATM) are reported descriptively only — p53
  targets such as CDKN1A also respond to p53-independent inputs, so a mean shift is not evidence
  against specificity.
- **H2 (exploratory, n = 2).** input_up KDs (NF1, TSC1) raise M1 without raising Z2.
- **M3 and the break-only set (exploratory).** Reported, no verdict.

## 8. Decision rule

- **Positive:** G2 ∧ H1 ∧ H1-spec, and N1 not violated → a loop-integrity readout exists in
  scRNA; proceed to build level 2 of the tool and replicate in a second screen (e.g. RPE1 or a
  stimulation dataset for input_up).
- **Null:** anything else → documented negative. The tool stays at level 1 (activity +
  readiness + audit), and the noise readout is reported as not detectable at Perturb-seq depth.
- Both outcomes are reported with every number above; no threshold changes after reading data.
  Deviations go to `docs/PREREG_loop_integrity_v1_DEVIATIONS.md`.

## 9. Settled before locking

1. Thresholds confirmed: 50 cells, fold_expr <= 0.50, D = 10th percentile, alpha 0.05.
2. K562 TP53 status: cited in §3, verified at run time.
3. Full download (disk available). The script reads the gem group from `obs['gem_group']` or,
   failing that, `obs['batch']`, and the knockdown from `obs['gene']` or `obs['perturbation']`;
   if neither exists it stops (that is a deviation, not a silent fallback).
4. UPR-ATF6 enters as a sensitivity analysis with fold_expr <= 0.60 (HSPA5 at 0.57). No verdict.

## 10. Validation of the locked code on synthetic data (before lock)

`tests/test_loop_integrity_prereg.py` plants known truth in simulated counts with the real
analysis set:

- **mean-only world** (break KDs raise the mean, no extra dispersion): decision NULL in 5/5
  seeds — the metric does not mistake a mean shift for a broken loop.
- **dispersion world** (break KDs raise dispersion 6x and the mean 1.4x): POSITIVE in 4/5
  seeds; the one miss was an N1 false alarm.
- **null calibration** (p53 KDs, no planted effect, 12 seeds, 60 units): Z2 sd 0.85, mean 0.08,
  Z2 >= 2.128 in 1.7 % per KD (about 5 % over three), |Z2| > 1.96 in 5 %. G2 pseudo-perturbation
  rates 1.4–8.3 % against the 10 % limit.

One design error was found and fixed during this validation, before lock: the draft N1 also
required the M1 CI to include 0 for all five p53 KDs; without multiplicity control that fails in
~23 % of null worlds, and a mean shift is not what the hypothesis forbids. N1 now tests Z2 of the
three break KDs only, Bonferroni-corrected.

## 11. What the run must return

`results/loop_integrity_v1/verdict.json`, `units.json`, `units.tsv` (unedited script output),
the figshare md5 and sha256 of the h5ad, the Cellosaurus TP53 annotation, software versions,
and the perceptome commit hash. Any deviation goes to
`docs/PREREG_loop_integrity_v1_DEVIATIONS.md` with the reason.
