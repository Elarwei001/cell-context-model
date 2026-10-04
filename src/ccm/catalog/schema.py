"""Context schema and identity normalisation.

A *context* is the combination of conditions that defines a cell's state: species, biological
identity (cell type or cell line), state (tissue, developmental stage, disease, donor, time point)
and intervention. Source adapters emit one row per context they contain, already aggregated
within that source; :mod:`ccm.catalog.build` merges rows across sources.
"""

from __future__ import annotations

import re

import pandas as pd

#: Columns every source adapter must return (extra columns are dropped).
CONTEXT_COLUMNS: list[str] = [
    "source",  # dataset or collection identifier, e.g. "census:2025-11-08:<dataset_id>"
    "species",
    "identity_kind",  # "cell_type" | "cell_line"
    "identity",  # canonical identity key (CL label for cell types; normalised name for lines)
    "cell_line_id",  # Cellosaurus accession where known, else ""
    "tissue",
    "tissue_general",
    "cell_type",
    "development_stage",
    "disease",
    "donor_id",  # source-scoped, never linked across sources
    "time_point",
    "intervention_kind",  # "none" | "genetic" | "chemical" | "ligand" | "other"
    "intervention_target",  # gene symbol, compound or ligand; "" when intervention_kind == "none"
    "intervention_detail",  # modality, dose, etc.
    "assay",
    "n_cells",  # single-cell count (0 for bulk sources)
    "n_samples",  # bulk profiles or signatures (0 for single-cell sources)
    "has_spatial",
    "has_kinetics",
    "licence",
]

#: Axes that define a context (the context key).
KEY_COLUMNS: list[str] = [
    "species",
    "identity_kind",
    "identity",
    "tissue",
    "development_stage",
    "disease",
    "donor_id",
    "time_point",
    "intervention_kind",
    "intervention_target",
    "intervention_detail",
]

INTERVENTION_KINDS = ("none", "genetic", "chemical", "ligand", "other")

#: Known cell lines with Cellosaurus accessions, keyed by normalised name.
CELLOSAURUS: dict[str, str] = {
    "K562": "CVCL_0004",
    "RPE1": "CVCL_4388",  # hTERT RPE-1
    "JURKAT": "CVCL_0065",
    "HEPG2": "CVCL_0027",
    "HCT116": "CVCL_0291",
    "HEK293T": "CVCL_0063",
    "A549": "CVCL_0023",
    "MCF7": "CVCL_0031",
    "A375": "CVCL_0132",
    "PC3": "CVCL_0035",
    "HELA": "CVCL_0030",
    "HT29": "CVCL_0320",
}

_ALIASES = {"HTERTRPE1": "RPE1", "RPE1HTERT": "RPE1", "HEK293FT": "HEK293T"}


def normalise_cell_line(name: str) -> str:
    """Canonical key for a cell-line name: upper case, alphanumerics only, known aliases folded.

    >>> [normalise_cell_line(n) for n in ("Hep G2", "HCT-116", "hTERT RPE-1")]
    ['HEPG2', 'HCT116', 'RPE1']
    """
    key = re.sub(r"[^A-Z0-9]", "", str(name).upper())
    return _ALIASES.get(key, key)


def empty_frame() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype=_dtype(c)) for c in CONTEXT_COLUMNS})


def _dtype(col: str) -> str:
    if col in ("n_cells", "n_samples"):
        return "int64"
    if col in ("has_spatial", "has_kinetics"):
        return "bool"
    return "object"


def conform(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` restricted and coerced to the context schema, filling missing columns."""
    out = pd.DataFrame(index=df.index)
    for col in CONTEXT_COLUMNS:
        if col in df.columns:
            out[col] = df[col]
        elif col in ("n_cells", "n_samples"):
            out[col] = 0
        elif col in ("has_spatial", "has_kinetics"):
            out[col] = False
        else:
            out[col] = ""
    for col in CONTEXT_COLUMNS:
        if _dtype(col) == "object":
            out[col] = out[col].fillna("").astype(str)
        elif _dtype(col) == "int64":
            out[col] = out[col].fillna(0).astype("int64")
        else:
            out[col] = out[col].fillna(False).astype(bool)
    out.loc[out["intervention_kind"] == "", "intervention_kind"] = "none"
    bad = set(out["intervention_kind"]) - set(INTERVENTION_KINDS)
    if bad:
        raise ValueError(f"unknown intervention_kind values: {sorted(bad)}")
    return out.reset_index(drop=True)
