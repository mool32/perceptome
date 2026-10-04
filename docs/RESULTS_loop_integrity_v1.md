# Results: loop integrity v1 (pre-registered)

**Decision (pre-registered rule): NULL.** A noise signature that distinguishes a broken feedback
loop from a removed input was **not** detected in K562 genome-wide Perturb-seq.

Pre-registration: `docs/PREREG_loop_integrity_v1.md`, locked at commit `e120d36`, tag
`prereg-loop-v1`. Run: 2026-10-04, 38 min, script unmodified (sha256 re-verified after the run).
Data: `K562_gwps_raw_singlecell_01.h5ad`, figshare md5 `887e3e6a…9546` matched, sha256
`b697ef7f…de2c`. Raw outputs: `results/loop_integrity_v1/` (unedited). Run report from the
executing agent: `results/loop_integrity_v1/REPORT_run.md`.

## 1. Verdict, as computed

| step | result |
|---|---|
| depth D | 5,223 UMI (6,165 cells below D dropped) |
| **G2** — metric vs sampling noise | 2.5 % of 199 pseudo-perturbations with \|Z2\| > 1.96 (limit 10 %) → **pass** |
| **G1** — readout reaches the loop | pass: HIF, JAK-STAT, SREBP · fail: ERK/MAPK, UPR-IRE1, UPR-PERK → 3 loops, exactly the minimum |
| **H1** — break > 0 and break > input_down | mean Z2(break) = 0.50 (n = 7), sign-flip p = 0.024 ✓; Mann–Whitney break vs input_down p = 0.24 ✗ → **fail** |
| **H1-spec** — break coefficient after \|M1\| | 0.64, p = 0.031 → pass |
| **N1** — p53 negative control (TP53-null K562) | no break KD with Z2 ≥ 2.128 → **not violated** |
| **Decision** | **NULL: decision rule not met** (H1 failed) |

Z2 values entering H1 — break: VHL −0.06, EGLN1 1.22, SOCS1 0.01, CISH 0.44, PTPN2 0.19,
INSIG1 0.74, INSIG2 0.93. input_down: HIF1A 0.58, ARNT −1.14, STAT5A 0.34, STAT5B 0.43, JAK2
−0.73, JAK1 −0.59, SCAP 0.92, SREBF2 0.19, MBTPS1 0.40, MBTPS2 1.20. Break knockdowns sit
slightly above zero, but so do the input knockdowns; the two classes are not separated.

TP53 status confirmed (Cellosaurus CVCL_0004: TP53 p.Gln136fs*13, homozygous), so N1 enters the
decision as planned.

## 2. Code–prose note (no effect on the decision)

The prose of §6 defines G1 on M1 and its CI. In the locked code a unit's M1 CI is only computed
when M2 is computable (≥ 3 output genes with mean ≥ 0.1), and G1 only counts such units. Effect:

- **ERK/MAPK** — input_down M1 = −0.08, −0.01, +0.12, +0.12 (median +0.06): fails G1 either way.
- **UPR-IRE1** — median input_down M1 −0.03, CI not available. Even if it had passed, its only
  break KD (DNAJB9) has no computable M2 and cannot enter H1.
- **UPR-PERK** — G1 saw only EIF2AK3 (M1 +0.06); ATF4 (M1 −0.54) had no M2. Had PERK passed G1,
  H1 would add PPP1R15A (Z2 −0.12) to break and EIF2AK3 (+0.10) to input_down, lowering the
  break mean. H1 still fails.

The decision is invariant to this ambiguity. Recorded here rather than in the deviations file,
because the code — which governs — was executed exactly.

## 3. Why the test had little to work with (descriptive, not part of the decision)

- **Readout panels were thin.** After thinning, 3–5 output genes per module passed the
  mean ≥ 0.1 filter; ERK and UPR-IRE1 had only 1–2, so M2 could not be computed. Per-unit M2
  bootstrap intervals span roughly ±2–3 log units, so effects of the size seen (Z2 ≈ 0.5)
  could not be resolved with 7 vs 10 units.
- **Many break knockdowns did not move the mean.** VHL (89 % knockdown) M1 = +0.05,
  EGLN1 −0.05, SOCS1 +0.03, CISH +0.10, PTPN2 −0.01, SPRY2/4 and SPRED1 ≈ +0.03. Either the
  loop is buffered (paralogs: MAPK3 for MAPK1, MAP2K2 for MAP2K1, SOCS family) or these
  outputs are saturated in K562. For those units the "break" class did not break anything
  measurable.

## 4. Exploratory observations (not pre-registered as hypotheses; reported for completeness)

- **Where a break clearly worked, dispersion did not consistently rise.** Break knockdowns
  with a clear mean increase: KEAP1 (M1 +1.00), HSPA5 (+1.48), INSIG1 (+0.65), CUL3 (+0.35),
  HSP90AB1 (+0.34). Their Z2: −0.33, +1.32, +0.74, −0.20, +1.21. The cleanest broken loop in
  the screen — NRF2 after KEAP1 loss, output doubled — shows **no** excess dispersion. At this
  depth and panel size this argues against the noise-compression readout rather than for it.
- **Mean-level specificity of the loop gene sets holds (pre-registered as descriptive).**
  Shared node MBTPS1/2 lowers both SREBP (M1 −0.56 / −0.88) and UPR-ATF6 outputs (−0.15 /
  −0.37); SCAP lowers SREBP (−0.50) but not UPR-ATF6 (+0.03); ATF6 lowers UPR-ATF6 (−0.13)
  but not SREBP (+0.04). Strong known breaks give strong mean shifts (KEAP1 → NRF2 targets,
  HSPA5 → ER chaperones). The level-1 (mean) readout of the loop layer is specific.
- **M3 coupling** (target–feedback correlation change) is near zero everywhere (|Δρ| ≤ 0.13
  outside one negative-control unit).

## 5. What this means for the tool

- Level 2 as specified — a noise signature of loop integrity readable from scRNA — **is not
  supported** by this pre-registered test. Per §8 of the pre-registration, the tool stays at
  level 1: activity + readiness + audit, with the loop layer used for **specific mean-level
  readouts** (break / input perturbation classes, shared-node specificity), which this run
  supports descriptively.
- A rescue on the same data (bigger regulons, different dispersion estimator, other thresholds)
  would be a new, post-hoc analysis of data already seen and must not be presented as a
  replication of this test. Any second attempt needs a new pre-registration and, preferably, an
  independent dataset with deeper per-gene counts.
