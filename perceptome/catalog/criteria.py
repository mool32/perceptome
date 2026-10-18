"""The four operational criteria of each module (paper Table S9b).

A module is defined by a signal-sensing component, a transcription-factor output with documented
targets, a negative-feedback element enabling cyclical operation, and receptor-field specificity.
The paper is explicit about what the table is: a record of how the 44 modules were *curated*, not a
verified property of the 206 genes as scored (Note S-N iv: the feedback criterion could not be
validated empirically). Its limits, all visible in this table:

* a feedback element is named for 38 of 44 modules; six are marked unclear (AR, PR, PPARgamma,
  NPAS4, PI3K/PTEN and Hippo, whose recorded entry was withdrawn);
* the receptor field is documented for 28 modules and blank for the 16 whose annotation was
  recorded at catalog construction;
* for 10 modules the named feedback element is not among the scored core genes (11 in the paper,
  which counts Hippo's withdrawn entry), so their score cannot read the feedback arm.

Gene lists and free-text cells are abbreviated in the print; the complete table is
``tableS9b_module_criteria.csv`` in the Zenodo record.
"""

import json
from functools import lru_cache
from pathlib import Path

import pandas as pd

from .modules import load_catalog

_FILE = Path(__file__).parent / "data" / "module_criteria_s9b.json"


@lru_cache(maxsize=1)
def _raw():
    with open(_FILE, encoding="utf-8") as f:
        return json.load(f)


def module_criteria():
    """Table S9b as a DataFrame indexed by module.

    Columns: ``sensor`` and ``feedback`` (lists as printed), ``tf_output`` (text as printed),
    ``tf_output_truncated``, ``feedback_named``, ``receptor_field_recorded``, ``confidence``
    (``rec.`` = recorded at construction, else high / medium / low), ``annotation``, and
    ``feedback_in_scored_genes`` (True / False for modules that name a feedback element, NaN otherwise).
    """
    cat = load_catalog()["modules"]
    rows = {}
    for m, v in _raw()["modules"].items():
        row = dict(v)
        core = set(cat[m]["core_genes"]) if m in cat else set()
        row["feedback_in_scored_genes"] = bool(set(v["feedback"]) & core) if v["feedback"] else float("nan")
        rows[m] = row
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("module")
