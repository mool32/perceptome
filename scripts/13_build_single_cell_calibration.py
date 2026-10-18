"""Build the single_cell_scaled reference calibration from the HPA single-cell-type table.

The paper projects per-cell data against a reference "re-derived under the single_cell_scaled
calibration — the same normalize_total(1e4)+log1p preprocessing applied to the tumor data — so
per-cell z-scores are commensurate" (Note S1); the pseudobulk calibration collapses single cells
onto low-complexity targets (Notes S-B, S-H). The calibration bundled with the package was built with
this script from the HPA single-cell-type table, which is not redistributed itself; rebuild it (or build
one from another HPA release) from your copy:

    # https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip  (unzip it first)
    python scripts/13_build_single_cell_calibration.py rna_single_cell_type.tsv

Steps
  1. Read the table; per cell type keep the 206 catalog genes and the total over all genes.
  2. Rebuild the pseudobulk reference (mean log1p(nCPM)) and compare it with the bundled 154 x 44
     matrix — an independent check that the table is read correctly and is the same release.
  3. Build the single_cell_scaled reference: each cell type scaled to 1e4 total, log1p, module means.
  4. Build the calibration (module mean / sample SD of that matrix, fixed 44 x 9 loadings, centroids).
     The sample SD (ddof=1) is the convention under which Table S7 is reproduced to the last printed
     digit; it differs from the population SD by one factor for every module, which changes no cosine.
  5. Check it against the paper: the z-scores of the megakaryocyte and cytotrophoblast over the 44
     modules must reproduce Table S7 (all modules but cGAS-STING, which the paper's own Note S-Q
     flags; see ``S7_KNOWN_DEVIATIONS``). Only then is the file written under the name that
     ``pct.load_calibration('single_cell_scaled')`` loads; otherwise it is written next to it as
     ``*.UNVERIFIED.json`` and the script exits with status 2.

Usage:  python scripts/13_build_single_cell_calibration.py <rna_single_cell_type.tsv> [--out PATH]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import perceptome as pct  # noqa: E402
from perceptome.convergence import check_calibration_against_table_s7  # noqa: E402
from perceptome.eigenspace.calibration import _loadings_frame, _SINGLE_CELL_FILE  # noqa: E402
from perceptome.eigenspace.project import _load_reference  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("hpa_tsv", help="HPA rna_single_cell_type.tsv")
    ap.add_argument("--out", default=str(_SINGLE_CELL_FILE), help="output calibration JSON")
    args = ap.parse_args(argv)

    cat = pct.load_catalog()
    genes = sorted({g for m in cat["modules"].values() for g in m["core_genes"]})
    print(f"Reading {args.hpa_tsv} ({len(genes)} catalog genes) ...")
    expr, library = pct.read_hpa_single_cell_type(args.hpa_tsv, genes=genes)
    print(f"  {expr.shape[1]} cell types, library sizes {library.min():.3g} .. {library.max():.3g}")

    bundled = pct.load_hpa_perceptivity()["R"].astype(float)
    missing = sorted(set(bundled.index) - set(expr.columns))
    if missing:
        sys.exit(f"The table lacks {len(missing)} of the 154 reference cell types, e.g. {missing[:5]}")
    expr, library = expr[bundled.index], library[bundled.index]

    pseudo = pct.hpa_reference_scores(expr, scaling="pseudobulk")[bundled.columns]
    diff = np.abs(pseudo.values - bundled.values)
    r = np.corrcoef(pseudo.values.ravel(), bundled.values.ravel())[0, 1]
    print(f"Pseudobulk rebuild vs bundled 154 x 44 matrix: max |diff| {np.nanmax(diff):.2e}, r = {r:.5f}"
          + ("" if np.nanmax(diff) < 1e-4 else "   <- not identical: a different HPA release?"))

    scaled = pct.hpa_reference_scores(expr, scaling="single_cell_scaled", library_size=library)[bundled.columns]
    ref = _load_reference()
    cal = pct.Calibration.from_scores(scaled, _loadings_frame(ref, list(ref["module_order"])),
                                      "single_cell_scaled", ddof=1)
    chk = check_calibration_against_table_s7(cal)
    print("Table S7 check: megakaryocyte r = {megakaryocyte_r:.4f}, cytotrophoblast r = {cytotrophoblast_r:.4f}; "
          "{n_exact} of {n_checked} printed values reproduced exactly once rounded, "
          "max |diff| over them = {max_abs_diff_checked:.4f}".format(**chk))
    for module, (dm, dc) in chk["known_deviations"].items():
        print(f"  not tested ({module}, flagged by Note S-Q): differs by {dm:+.3f} / {dc:+.3f}")

    out = Path(args.out)
    if chk["matches"]:
        cal.save(out)
        print(f"Matches the paper. Wrote {out}")
        return 0
    unverified = out.with_name(out.stem + ".UNVERIFIED.json")
    cal.save(unverified)
    print(f"Does NOT reproduce Table S7 (needs max |diff| <= 0.006 outside the known deviations). Wrote {unverified}; "
          "pct.load_calibration('single_cell_scaled') will not load it.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
