"""Re-express the bundled reference shift vectors in the standardized coordinates.

Versions <= 0.2.2 projected module deltas as ``delta @ loadings``. The paper's coordinate
system is standardized (module / reference SD, then loadings), and ``pct.project`` /
``pct.project_shift`` now implement it, so the eigenspace-side vectors stored with the
bundled references have to live in the same convention or ``compare_to_references``
would compare a standardized cell position with an unscaled direction.

Every reference stores its 44-module deltas, which makes the conversion exact and
reproducible without the raw scRNA-seq data:

  attractor_v1.json      attractor_direction_eigenspace     <- attractor_direction_modules
  disease_vectors.json   per_cell_type[*].eigenspace_vector,
                         mean_vector_12D                    <- module_deltas / module_delta_mean_across_groups
  aging_reference.json   tissue_vectors[*].vector_PCs,
                         inflammaging_direction, collapse_direction  <- tissue_vectors[*].module_deltas

Idempotent: running it again leaves the files unchanged. The legacy (unscaled) vector is
``delta @ loadings`` and can always be recovered from the stored module deltas.

Usage:  python scripts/12_restandardize_reference_vectors.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import perceptome as pct  # noqa: E402

DATA = Path(__file__).resolve().parents[1] / "perceptome" / "reference" / "data"
CONVENTION = "standardized"
NOTE = ("Eigenspace vectors are module deltas divided by the reference SD and projected through the "
        "fixed loadings (pct.project_shift), i.e. the paper's standardized coordinates. Versions "
        "<= 0.2.2 stored the unscaled projection (delta @ loadings).")


def _vec(module_deltas):
    return pct.project_shift(pd.Series(module_deltas, dtype=float)).values


_STYLE = {}


def _read(name):
    """Load a JSON file and remember how it was written, so the diff stays minimal."""
    raw = (DATA / name).read_text(encoding="utf-8")
    _STYLE[name] = {"ensure_ascii": "\\u" in raw, "newline": raw.endswith("\n")}
    return json.loads(raw)


def _write(name, obj):
    style = _STYLE[name]
    text = json.dumps(obj, indent=2, ensure_ascii=style["ensure_ascii"])
    (DATA / name).write_text(text + ("\n" if style["newline"] else ""), encoding="utf-8")
    print(f"  wrote {name}")


def attractor():
    d = _read("attractor_v1.json")
    d["attractor_direction_eigenspace"] = _vec(d["attractor_direction_modules"]).tolist()
    d["eigenspace_convention"] = CONVENTION
    d["eigenspace_convention_note"] = NOTE
    _write("attractor_v1.json", d)


def disease():
    d = _read("disease_vectors.json")
    for name, rec in d.items():
        if name.startswith("_"):
            continue
        for group in rec["per_cell_type"].values():
            v = _vec(group["module_deltas"])
            group["eigenspace_vector"] = v.tolist()
            group["magnitude"] = float(np.linalg.norm(v))
        rec["mean_vector_12D"] = _vec(rec["module_delta_mean_across_groups"]).tolist()
    d["_meta"]["eigenspace_convention"] = CONVENTION
    d["_meta"]["eigenspace_convention_note"] = NOTE
    _write("disease_vectors.json", d)


def aging():
    d = _read("aging_reference.json")
    for tv in d["tissue_vectors"].values():
        v = _vec(tv["module_deltas"])
        tv["vector_PCs"] = v.tolist()
        tv["magnitude"] = float(np.linalg.norm(v))
        tv["pc1_direction"] = "up" if v[0] > 0 else "down"
    d["inflammaging_direction"] = d["tissue_vectors"]["blood"]["vector_PCs"]
    d["collapse_direction"] = d["tissue_vectors"]["bone_marrow"]["vector_PCs"]
    d["eigenspace_convention"] = CONVENTION
    d["eigenspace_convention_note"] = NOTE
    _write("aging_reference.json", d)


def main():
    print("Re-expressing bundled reference vectors in standardized coordinates")
    attractor()
    disease()
    aging()


if __name__ == "__main__":
    main()
