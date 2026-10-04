"""LINCS L1000 (CMap LINCS 2020 release): bulk signatures across hundreds of cell lines.

Uses ``siginfo_beta.txt`` (one row per consensus signature) and ``cellinfo_beta.txt`` (Cellosaurus
accessions). Counts are signatures and replicate profiles (``n_samples``), not single cells.
Perturbation types map to intervention kinds:
trt_sh / trt_sh.cgs / trt_xpr / trt_oe / trt_si -> genetic; trt_cp -> chemical;
trt_lig -> ligand; ctl_* -> none. trt_sh.css (shRNA seed-matched controls) and trt_aby
(antibodies) are kept as "other".
"""

from __future__ import annotations

import pandas as pd

from ..schema import normalise_cell_line

BASE = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020"
LICENCE = "LINCS data release policy (open)"
GENETIC = {"trt_sh": "shRNA", "trt_sh.cgs": "shRNA (consensus)", "trt_xpr": "CRISPR KO",
           "trt_oe": "overexpression", "trt_oe.mut": "overexpression (mutant)",
           "trt_si": "siRNA"}


def _kind(pert_type: str) -> str:
    if pert_type in GENETIC:
        return "genetic"
    if pert_type == "trt_cp":
        return "chemical"
    if pert_type == "trt_lig":
        return "ligand"
    if pert_type.startswith("ctl_"):
        return "none"
    return "other"


def l1000_contexts(siginfo_path: str, cellinfo_path: str) -> pd.DataFrame:
    sig = pd.read_csv(
        siginfo_path,
        sep="\t",
        usecols=["pert_type", "cmap_name", "cell_iname", "pert_idose", "pert_itime", "nsample"],
        dtype=str,
        low_memory=False,
    ).fillna("")
    cells = pd.read_csv(
        cellinfo_path, sep="\t", usecols=["cell_iname", "cellosaurus_id"], dtype=str
    )
    cellosaurus = dict(zip(cells["cell_iname"], cells["cellosaurus_id"].fillna("")))
    sig["intervention_kind"] = sig["pert_type"].map(_kind)
    ctrl = sig["intervention_kind"] == "none"
    sig["intervention_target"] = sig["cmap_name"].where(~ctrl, "")
    detail = sig["pert_type"].map(GENETIC).fillna("")
    dose = sig["pert_idose"].str.strip()
    sig["intervention_detail"] = (detail + " " + dose).str.strip().where(~ctrl, "")
    sig["time_point"] = sig["pert_itime"].str.strip()
    sig["n_samples"] = pd.to_numeric(sig["nsample"], errors="coerce").fillna(1).astype("int64")
    keys = ["cell_iname", "intervention_kind", "intervention_target", "intervention_detail",
            "time_point"]
    df = sig.groupby(keys, as_index=False)["n_samples"].sum()
    df["identity"] = df["cell_iname"].map(normalise_cell_line)
    df["cell_line_id"] = df["cell_iname"].map(cellosaurus).fillna("")
    df["source"] = "lincs-l1000:2020"
    df["species"] = "human"
    df["identity_kind"] = "cell_line"
    df["assay"] = "L1000 (bulk, ~978 landmark genes)"
    df["licence"] = LICENCE
    return df
