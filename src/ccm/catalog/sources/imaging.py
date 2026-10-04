"""Image-based perturbation screens (Cell Painting and related multiplexed stains).

Counts are wells (``n_samples``); per-well cell counts are not in the public metadata.

- JUMP Cell Painting (cpg0016, Cell Painting Gallery): U2OS cells; compounds, CRISPR knockouts
  and ORF overexpression. Metadata from github.com/jump-cellpainting/datasets (metadata/*.csv.gz).
- RxRx3-core (Recursion): HUVEC (human umbilical vein endothelial cells); CRISPR knockouts and
  compounds at several concentrations. Metadata from Hugging Face recursionpharma/rxrx3-core.
"""

from __future__ import annotations

import os

import pandas as pd

JUMP_METADATA = "https://github.com/jump-cellpainting/datasets/raw/main/metadata"
JUMP_LICENCE = "CC0-1.0"
JUMP_DMSO = "JCP2022_033924"
JUMP_CRISPR_CONTROLS = {"no-guide", "non-targeting"}

RXRX3_METADATA = (
    "https://huggingface.co/datasets/recursionpharma/rxrx3-core/resolve/main/"
    "metadata_rxrx3_core.csv"
)
RXRX3_LICENCE = "Recursion RxRx3-core licence (restrictive; see dataset)"
HUVEC_CELL_TYPE = "endothelial cell of umbilical vein"  # Cell Ontology CL:0002618


def _read(path_or_url: str, **kw) -> pd.DataFrame:
    return pd.read_csv(path_or_url, low_memory=False, **kw)


def jump_contexts(metadata_dir: str | None = None) -> pd.DataFrame:
    """JUMP cpg0016: one row per (perturbation) in U2OS, counting wells."""

    def src(name: str) -> str:
        local = os.path.join(metadata_dir, f"{name}.csv.gz") if metadata_dir else ""
        return local if local and os.path.exists(local) else f"{JUMP_METADATA}/{name}.csv.gz"

    wells = _read(src("well"))
    compound = _read(src("compound"), usecols=["Metadata_JCP2022", "Metadata_InChIKey"])
    crispr = _read(src("crispr"), usecols=["Metadata_JCP2022", "Metadata_Symbol"])
    orf = _read(src("orf"), usecols=["Metadata_JCP2022", "Metadata_Symbol", "Metadata_pert_type"])

    kind, target, detail = {}, {}, {}
    for jcp, key in zip(compound["Metadata_JCP2022"], compound["Metadata_InChIKey"]):
        ctrl = jcp == JUMP_DMSO
        kind[jcp] = "none" if ctrl else "chemical"
        target[jcp] = "" if ctrl else str(key)
        detail[jcp] = "" if ctrl else "compound (InChIKey)"
    for jcp, sym in zip(crispr["Metadata_JCP2022"], crispr["Metadata_Symbol"]):
        ctrl = str(sym) in JUMP_CRISPR_CONTROLS
        kind[jcp] = "none" if ctrl else "genetic"
        target[jcp] = "" if ctrl else str(sym)
        detail[jcp] = "" if ctrl else "CRISPR KO"
    for jcp, sym, pt in zip(orf["Metadata_JCP2022"], orf["Metadata_Symbol"],
                            orf["Metadata_pert_type"]):
        ctrl = str(pt) != "trt"
        kind[jcp] = "none" if ctrl else "genetic"
        target[jcp] = "" if ctrl else str(sym)
        detail[jcp] = "" if ctrl else "ORF overexpression"

    w = wells[wells["Metadata_JCP2022"].isin(kind)]
    counts = w.groupby("Metadata_JCP2022").size()
    df = pd.DataFrame({"jcp": counts.index, "n_samples": counts.to_numpy()})
    df["intervention_kind"] = df["jcp"].map(kind)
    df["intervention_target"] = df["jcp"].map(target)
    df["intervention_detail"] = df["jcp"].map(detail)
    df = df.groupby(["intervention_kind", "intervention_target", "intervention_detail"],
                    as_index=False)["n_samples"].sum()
    df["source"] = "jump-cpg0016"
    df["modality"] = "imaging"
    df["species"] = "human"
    df["identity_kind"] = "cell_line"
    df["identity"] = "U2OS"
    df["cell_line_id"] = "CVCL_0042"
    df["assay"] = "Cell Painting"
    df["licence"] = JUMP_LICENCE
    return df


def rxrx3_contexts(metadata_path: str | None = None) -> pd.DataFrame:
    """RxRx3-core: one row per (perturbation, concentration) in HUVEC, counting wells."""
    m = _read(metadata_path or RXRX3_METADATA, dtype={"gene": str, "SMILES": str})
    crispr = m["perturbation_type"] == "CRISPR"
    gene = m["gene"].fillna("").astype(str).str.strip()
    is_ctrl = crispr & (
        gene.eq("") | gene.str.upper().isin({"EMPTY_CONTROL", "INTRON", "CONTROL"})
        | m["well_type_label"].astype(str).str.contains("Intron controls")
    )
    m["intervention_kind"] = "chemical"
    m.loc[crispr, "intervention_kind"] = "genetic"
    m.loc[is_ctrl, "intervention_kind"] = "none"
    m["intervention_target"] = m["treatment"].astype(str).where(~crispr, gene)
    m.loc[is_ctrl, "intervention_target"] = ""
    conc = m["concentration"].map(lambda c: "" if pd.isna(c) else f"{c:g} uM")
    m["intervention_detail"] = ("compound " + conc).where(~crispr, "CRISPR KO")
    m.loc[is_ctrl, "intervention_detail"] = ""
    df = (
        m.groupby(["intervention_kind", "intervention_target", "intervention_detail"])
        .size()
        .rename("n_samples")
        .reset_index()
    )
    df["source"] = "rxrx3-core"
    df["modality"] = "imaging"
    df["species"] = "human"
    df["identity_kind"] = "cell_type"
    df["identity"] = HUVEC_CELL_TYPE
    df["cell_type"] = HUVEC_CELL_TYPE
    df["tissue"] = "umbilical vein"
    df["assay"] = "multiplexed fluorescence imaging (Cell Painting-like)"
    df["licence"] = RXRX3_LICENCE
    return df
