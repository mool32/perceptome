# Exploratory notes: dose-response, buffering, RPE1 inventory

**Status: exploratory.** K562_gwps data had already been used by the locked loop-integrity
test; nothing here is confirmatory. Two definitions were tried for the locality contrast
(§3), so its p-values are not corrected for that choice. Outputs: `results/explore_dose_buffering/`
(parts A, B from commit b9deef3; part C rerun from 54edab3 after dropping 73 genes with +inf
in Replogle's normalized pseudobulk), `results/rpe1_candidates/`.

## 1. Per-cell dose-response (parts A, B)

- Only 13/112 knocked-down genes are readable per cell (control mean >= 1 UMI, >= 50 % nonzero).
- Among readable cases with a defined role, the sign of the within-population correlation
  (residual knocked-down gene vs module output) matched the role in 10/11 cases: inhibitors /
  feedback nodes negative (HSPA8 −0.15, HSP90AB1, HSP90AA1, HSPA5, CUL3), inputs positive
  (GRB2, MAPK1, XBP1, ATF4, RHEB +0.22). Exception FKBP5 (GR inactive without ligand).
- **Second look: this is close to a tautology.** Under a monotone dose-response the
  within-population sign equals the sign of the knockdown's mean effect, which is already
  known from the mean. The only non-trivial part — feedback genes correlate *positively* with
  output in unperturbed cells (co-induction) but act negatively — is a known property of
  negative feedback. Not pursued as a new readout.
- **RPE1 cannot test it anyway:** 3 eligible units, all of the break class (CUL3, HSPA5 ×2),
  no inputs. RPE1 data remain unseen beyond control-cell levels of knocked-down genes.

## 2. Evolutionary cohorts (Paper 7) — no support visible

Response locality (own-module activity genes vs whole transcriptome, median over efficient
knockdowns): ancient-unit HSF1 1.96, AMPK 0.91, ERK/MAPK 0.91, PI3K/PTEN 0.48; young-anti-unit
p53 1.00, cGAS-STING 0.33, NFAT 1.00 (1 KD), VDR (1 KD, no locality). No separation; the
cohort test is not feasible in K562 (inactive young modules, n = 1).

## 3. Unexpected pattern: response locality of closed homeostatic loops

Locality correlates with baseline activity of the module (Spearman 0.34, p = 0.04), so it was
residualised on log baseline. After that:

- curated tier A (clean transcriptional loops active in K562) vs the rest: residual median
  +0.46 vs −0.29, Mann–Whitney p = 0.003 (tier A partly encodes K562 activity — circular);
- modules with a curated transcriptional `feedback_core` vs without: +0.24 vs −0.42,
  p = 0.031 (definition independent of K562).

Highest: SREBP 7.4 (knockdown of sterol-pathway enzymes and of INSIG/SCAP stays within the
SREBP program), UPR-ATF6 3.5, NRF2 3.0, HIF 2.5, HSF1 2.0, UPR-PERK 1.9.

Candidate hypothesis (**homeostatic closure**): a perturbation inside a loop that senses its
own output is compensated within that loop, so the transcriptional response stays local.
Confounds to remove in any confirmatory test: detectability of highly expressed activity
genes; knockdown lists that include the module's own targets (partly definitional); one cell
line; two definitions tried here. Requires a new pre-registration on an independent
genome-scale Perturb-seq dataset.
