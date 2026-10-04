# Amendment 1 to pre-registration v2 (data access only)

Date: 2026-10-04. Made **before any knockdown or readout expression value of X-Atlas/Orion was read**.
State at amendment: `--inspect` on both lines (column names, category labels, X layout); one aborted
build of HCT116 that fetched raw counts of sampled control cells over the network for ~30 min and
crashed (connection reset) before computing or printing anything.

Nothing in the hypothesis, metric, parameters, decision rules or `scripts/18_...` changes.

## 1. `read_column` handles anndata nullable arrays (script 17)

The atlas stores some obs string columns in the anndata >= 0.12 nullable encoding (group with
`values` and `mask`, mask True = missing). Script 17 did not decode it. Added one branch that returns
the values with missing entries as "" (strings) or 0 (numbers). Categorical and plain columns are
decoded exactly as before. The operator had applied an equivalent local fix (reported as a deviation);
this commit replaces it.

## 2. Sequential subset extraction (new `scripts/17a_extract_subset.py`)

Files exist only on Figshare+ (HCT116 209 GB, HEK293T 350 GB; the HF repository holds per-batch parquet,
not the h5ad). Script 17 reads needed rows by random access; over HTTP each scattered row costs a full
cache block (~160 GB for the controls alone), so the build cannot finish. 17a makes one sequential,
resumable pass over the CSR matrix and keeps only the rows script 17 would read — the control sample
drawn with script 17's own rule and seed, and every filtered cell of the analysis knockdowns — into a
small local h5ad. Script 17 then runs on that file with the same flags. A test
(`test_extractor_subset_gives_identical_pseudobulk`) shows the pseudobulk from the subset equals the
pseudobulk from the full file, including control sampling and batch normalisation.
