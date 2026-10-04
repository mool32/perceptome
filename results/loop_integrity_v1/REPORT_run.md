# Loop integrity v1 — run report

## Lock verification
- commit checked out: `e120d361e8776ba1d286898d70be8b36e77eb135` (detached HEAD at `e120d36`, "Lock pre-registration v1: loop integrity (noise signature of a broken loop)")
- PREREG_LOCK sha256 check: **all OK (6/6)** — verified before the run and re-verified after it, unchanged: `docs/PREREG_loop_integrity_v1.md`, `scripts/14_loop_integrity_prereg_v1.py`, `data/evidence/replogle2022_K562_gwps_contrasts.tsv`, `perceptome/catalog/data/loops_v04.json`, `perceptome/catalog/loops.py`, `tests/test_loop_integrity_prereg.py`. (`sha256sum` is absent on macOS; `shasum -a 256 -c` was used, as the brief allows.)
- pytest: **87 passed** (40.4 s)
- tag prereg-loop-v1 pushed: **yes** — `* [new tag] prereg-loop-v1 -> prereg-loop-v1` to `https://github.com/mool32/perceptome`, annotated on `e120d36`

## Data
- file: **K562_gwps_raw_singlecell_01.h5ad** | 65,830,941,948 bytes (65.83 GB) | figshare md5 `887e3e6a8c8df6eadf7a3030a53c9546` (`supplied_md5`; `computed_md5` is empty for this file) | **md5 match: yes** | sha256 `b697ef7fedcec2972ec334608f86f2a29cb25e70c81b6725195f4cda60a9de2c`
- obs columns: `gem_group, gene, gene_id, transcript, gene_transcript, sgID_AB, mitopercent, UMI_count, z_gemgroup_UMI, core_scale_factor, core_adjusted_UMI_count` (obs index = `cell_barcode`); **knockdown column: `gene`**; **gem-group column: `gem_group`**
- non-targeting cells: **75,328** (label is exactly `non-targeting`, a single matching category)
- X: **dense**, dtype `float32`, shape **(1,989,578, 8,248)**; integer check on 1 control cell × 50 genes (obs row 95, the first non-targeting cell): **pass** — `np.all(x == np.round(x))` = True, range 0–17, values 0.0/1.0/2.0/… → raw UMI counts
- var gene symbols from: **`gene_name`** (symbols: `LINC01409, LINC01128, NOC2L, KLHL17, HES4, …`). The var index itself is Ensembl (`ENSG…`), but the locked script reads `a.var["gene_name"]` when present (script lines 97–98), so the §4.4 stop condition did not apply.

## TP53 (Cellosaurus CVCL_0004)
Verbatim, from `https://api.cellosaurus.org/cell-line/CVCL_0004?format=txt` (line 1517):

```
CC   Sequence variation: Mutation; HGNC; HGNC:11998; TP53; Simple; p.Gln136fs*13 (c.406_407insC); Zygosity=Homozygous (PubMed=17088437; PubMed=18277095).
```

Also present (line 1516), recorded for context:

```
CC   Sequence variation: Gene fusion; HGNC; HGNC:76; ABL1 + HGNC; HGNC:1014; BCR; Name(s)=BCR-ABL1, BCR-ABL; Note=BCR exon 14 fused to ABL1 exon 2 (b3a2 transcript) (PubMed=10071072; PubMed=10071072; PubMed=12506034; PubMed=20809971).
```

## Run
- command:
  `python scripts/14_loop_integrity_prereg_v1.py --h5ad ~/data/replogle2022/K562_gwps_raw_singlecell_01.h5ad --out results/loop_integrity_v1 2>&1 | tee results/loop_integrity_v1_run.log`
- start 2026-10-04 08:12:45 / end 2026-10-04 08:50:54, wall clock **38 min 09 s**, exit code 0
- restarts: **none**
- python 3.12.0 / numpy 2.5.3 / scipy 1.18.1 / anndata 0.13.4 / h5py 3.16.0

## Verdict

Contents of `results/loop_integrity_v1/verdict.json`, unedited:

```json
{
 "depth_D": 5223,
 "n_cells_dropped_below_D": 6165,
 "G2": {
  "frac_abs_Z2_gt_1.96": 0.02512562814070352,
  "n": 199,
  "pass": true
 },
 "G1": {
  "ERK/MAPK": {
   "median_M1_input_down": NaN,
   "any_ci_below_0": false,
   "pass": false
  },
  "HIF": {
   "median_M1_input_down": -0.22564039360127436,
   "any_ci_below_0": true,
   "pass": true
  },
  "JAK-STAT": {
   "median_M1_input_down": -0.08441566632873873,
   "any_ci_below_0": true,
   "pass": true
  },
  "SREBP": {
   "median_M1_input_down": -0.5290524146662746,
   "any_ci_below_0": true,
   "pass": true
  },
  "UPR-IRE1": {
   "median_M1_input_down": NaN,
   "any_ci_below_0": false,
   "pass": false
  },
  "UPR-PERK": {
   "median_M1_input_down": 0.05534736810615705,
   "any_ci_below_0": false,
   "pass": false
  }
 },
 "loops_passing_G1": [
  "HIF",
  "JAK-STAT",
  "SREBP"
 ],
 "H1": {
  "n_break": 7,
  "n_input_down": 10,
  "mean_Z2_break": 0.49636964107643916,
  "p_sign_flip": 0.023897610238976102,
  "p_mann_whitney": 0.23735088440970795,
  "pass": false
 },
 "H1_spec": {
  "break_coef": 0.6374850607789154,
  "p_one_sided": 0.03078230142052135,
  "pass": true
 },
 "N1": {
  "violated_by": [],
  "not_violated": true
 },
 "decision": "NULL: decision rule not met"
}
```

## Deviations

1. **Disk space had to be freed before the data could be downloaded.** The volume had 33 GB free against the 65.83 GB the file needs. On the user's explicit instruction and after his confirmation of an itemised list, five public re-downloadable dataset files (49.9 GB) were deleted from `research/perceptual_modules`: `cmap_data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx` (23.38 GB), `paper6/data/perturbseq/ReplogleWeissman2022_K562_gwps.h5ad` (8.81 GB), `paper1/oscilatory/data/tms/tms_droplet.h5ad` (8.24 GB), `paper1/oscilatory/data/tms/tms_facs.h5ad` (4.80 GB), `paper6/data/weinreb2020/GSE140802_RAW.tar` (4.67 GB). No results, no derived files, no paper-3 pipeline inputs were touched. This is outside the repository and does not affect the run.
2. **`shasum -a 256 -c` instead of `sha256sum`** for the lock check — the brief provides for this on macOS.
3. The download used `curl -C -` (resume) and was verified by exact byte size and then by md5 against figshare, as the brief's §3 allows.
4. Nothing in the repository was modified. `git status` is clean apart from untracked `.venv/` and `results/`. Nothing was committed; the only push was the tag in §2.3.

## Confirmation
- expression values were not read before the run except the 1×50 integer check: **yes** — exactly one slice, `h['X'][95, :50]`, on the first `non-targeting` cell; no other expression value was read at any point
- the script was not modified: **yes** — sha256 `416019e1964379120839eb2c1c0d5cfb7c8e8f2505eb9e5da56e24a3cd593997`, re-verified OK after the run
- no additional analyses or plots were made after the run: **yes** — the outputs are returned exactly as the script wrote them; nothing was recomputed, re-read differently, re-run with other parameters, or plotted

## Attached
`results/loop_integrity_v1/verdict.json`, `results/loop_integrity_v1/units.json`, `results/loop_integrity_v1/units.tsv`, `results/loop_integrity_v1_run.log` (64 lines — the script's whole stdout is the verdict JSON).
