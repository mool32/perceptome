"""Gene-symbol aliases for catalog genes whose HGNC symbol has changed.

Note S-Q of the paper audited the scored panel against 17 cohorts and found catalog genes
that matched neither the Human Protein Atlas reference nor most cohorts, for two reasons:

* deprecated symbols — MB21D1 and TMEM173 (cGAS-STING; now CGAS and STING1) and ARNTL
  (Circadian; now BMAL1). Matching by exact symbol therefore scored cGAS-STING on IRF3, IFNB1
  and CXCL10 only — without either of its sensors — and Circadian on four genes of five.
* a loader that stripped any trailing ``-<digits>`` (meant to undo ``var_names_make_unique``)
  also turned the real symbol NKX3-1 into NKX3, so AR was scored on four genes of five in every
  cohort. This package never strips suffixes; see ``perceptome.convergence.prepare_expression``.

The catalog keeps the symbols exactly as published (Table S9). ``resolve_aliases=True`` in
``score_modules`` additionally matches the current symbol, so a dataset annotated with either
spelling scores the same five genes. It is off by default because the bundled reference matrices
and eigenspace were built by exact matching on HPA (current symbols): cGAS-STING on three genes,
Circadian on four. Re-deriving the reference with the corrected catalog is the one step the
paper lists as not completed.
"""

# catalog symbol -> current HGNC symbol(s)
_PAPER_AUDITED = {
    "MB21D1": ("CGAS",),
    "TMEM173": ("STING1",),
    "ARNTL": ("BMAL1",),
}

# Renamed symbols in the *activity* gene lists (not part of the paper's audit of the scored panel).
_ACTIVITY_LISTS = {
    "CTGF": ("CCN2",),
    "CYR61": ("CCN1",),
    "ACPP": ("ACP3",),
    "GILZ": ("TSC22D3",),
}

ALIASES = {**_PAPER_AUDITED, **_ACTIVITY_LISTS}


def aliases_of(gene):
    """Alternative symbols under which a catalog gene may appear in a dataset (may be empty)."""
    return ALIASES.get(gene, ())


def resolve_genes(genes, available):
    """Match catalog genes against the gene names of a dataset, falling back to aliases.

    Parameters
    ----------
    genes : sequence of str
        Catalog gene symbols.
    available : set-like of str
        Gene names present in the dataset.

    Returns
    -------
    (found, missing, aliased)
        ``found``   names to read from the dataset (the catalog symbol, or its alias)
        ``missing`` catalog genes with no match
        ``aliased`` {catalog symbol: alias used}
    """
    found, missing, aliased = [], [], {}
    for g in genes:
        if g in available:
            if g not in found:
                found.append(g)
            continue
        alt = next((a for a in aliases_of(g) if a in available), None)
        if alt is None:
            missing.append(g)
        else:
            aliased[g] = alt
            if alt not in found:
                found.append(alt)
    return found, missing, aliased
