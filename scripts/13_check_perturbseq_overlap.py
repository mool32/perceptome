"""Which loops are testable in a Perturb-seq screen?

Reads only metadata (backed mode), never the expression matrix.

Usage:
  python scripts/13_check_perturbseq_overlap.py --h5ad K562_essential_raw_singlecell_01.h5ad
  python scripts/13_check_perturbseq_overlap.py --perts perturbed.txt --readout measured.txt

--h5ad    perturbations from obs['gene'] (Replogle convention) or obs['perturbation'];
          readout genes from var['gene_name'] or var_names; cells counted per knockdown.
--perts   one knocked-down gene symbol per line; --readout one measured gene per line.
"""

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import perceptome as pct  # noqa: E402


def load_h5ad(path):
    import anndata as ad
    a = ad.read_h5ad(path, backed="r")
    col = "gene" if "gene" in a.obs else "perturbation"
    cells = collections.Counter(a.obs[col].astype(str))
    readout = set(a.var["gene_name"].astype(str)) if "gene_name" in a.var else set(map(str, a.var_names))
    return cells, readout


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad")
    p.add_argument("--perts")
    p.add_argument("--readout")
    p.add_argument("--min_cells", type=int, default=50)
    p.add_argument("--tiers", default="ABCD")
    a = p.parse_args()

    if a.h5ad:
        cells, readout = load_h5ad(a.h5ad)
    elif a.perts:
        cells = collections.Counter({l.strip(): a.min_cells for l in open(a.perts) if l.strip()})
        readout = {l.strip() for l in open(a.readout)} if a.readout else set()
    else:
        p.error("give --h5ad or --perts")

    out = sys.stdout
    out.write("module\ttier\tbreak\tinput_down\tinput_up\tfeedback_measured\toutput_measured\ttestable\n")
    for name in pct.list_loops(tuple(a.tiers)):
        e = pct.get_loop(name)
        present = {k: [f"{g}({cells[g]})" for g in v if cells.get(g, 0) >= a.min_cells]
                   for k, v in e["perturbations"].items()}
        fb = [g for g in e["feedback_core"] if g in readout]
        output = pct.loop_output_genes(name)
        tg = [g for g in output if g in readout]
        has_contrast = bool(present["break"] or present["input_up"]) and bool(present["input_down"])
        enough_readout = len(tg) >= 3 and (fb or not e["feedback_readable"] or e["feedback_is_output"])
        testable = has_contrast and (enough_readout if readout else True)
        out.write("\t".join([
            name, e["tier"], ",".join(present["break"]), ",".join(present["input_down"]),
            ",".join(present["input_up"]), f"{len(fb)}/{len(e['feedback_core'])}",
            f"{len(tg)}/{len(output)}", "YES" if testable else "no",
        ]) + "\n")


if __name__ == "__main__":
    main()
