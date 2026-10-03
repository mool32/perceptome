# Evidence tables (metadata only)

Source: Replogle et al. 2022, figshare article 20029387 (v1), pseudobulk files
`K562_essential_normalized_bulk_01.h5ad`, `K562_gwps_normalized_bulk_01.h5ad`,
`rpe1_normalized_bulk_01.h5ad`. Extracted 2026-10 by reading only the `obs` and
`var` groups (expression matrix never read). sha256 of the source files:

- K562_essential_normalized_bulk_01.h5ad  c1ca6456c9c9f1aa2b02c496eb64d1dc3e6a852edbd744d682b8d2c95fd36829
- rpe1_normalized_bulk_01.h5ad            a3c5bfd0f15d63938bc80c9b8874b9cd761e3a23caf5ffe7966bae4e887ec89d
- K562_gwps_normalized_bulk_01.h5ad       37e48c474d8b5dead4151f96ea8f5fe7bbe6beb10eeea48685b740c3f74490a2

Files:
- `replogle2022_testability.tsv` — per loop and screen: perturbation contrast present (break or
  input_up, and input_down, each >= 50 cells), >= 3 loop targets measured, feedback measured
  (rule of `scripts/13_check_perturbseq_overlap.py` at commit 894d008).
- `replogle2022_K562_gwps_contrasts.tsv` — every perturbation gene of every loop in K562_gwps:
  presence, `num_cells_filtered` summed over guide groups, author `fold_expr` (residual
  fraction of the targeted transcript; copied verbatim, not recomputed).

Symbol mapping was done by Ensembl id against HGNC; no loop-layer gene needed an
alias. Three homonyms were rejected (KCNA4~HK1, PMEPA1~STAG1, WWTR1~TAZ).
Control cells in K562_gwps: 75,328 (matches the published count).
