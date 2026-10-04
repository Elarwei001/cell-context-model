"""Tahoe-100M (Tahoe Therapeutics): ~100M cells, 50 cancer cell lines x ~380 drugs x 3 doses.

Reads the per-cell ``obs_metadata.parquet`` from Hugging Face with column projection (cell line,
drug and concentration only) and aggregates cell counts per (cell line, drug, dose).
DMSO is treated as the vehicle control (``intervention_kind="none"``).
"""

from __future__ import annotations

import ast

import pandas as pd

from ..schema import normalise_cell_line

REPO = "datasets/tahoebio/Tahoe-100M"
LICENCE = "CC0-1.0"
VEHICLES = {"DMSO_TF", "DMSO"}


def _dose(drugconc: str) -> str:
    try:
        parsed = ast.literal_eval(drugconc)
        return ";".join(f"{c} {u}" for _, c, u in parsed)
    except (ValueError, SyntaxError):
        return ""


def tahoe_contexts(verbose: bool = True) -> pd.DataFrame:
    import pyarrow.parquet as pq
    from huggingface_hub import HfFileSystem

    fs = HfFileSystem()
    with fs.open(f"{REPO}/metadata/cell_line_metadata.parquet") as fh:
        lines = pq.read_table(fh, columns=["cell_name", "Cell_ID_Cellosaur"]).to_pandas()
    cellosaurus = dict(zip(lines["cell_name"], lines["Cell_ID_Cellosaur"]))
    counts: pd.Series | None = None
    with fs.open(f"{REPO}/metadata/obs_metadata.parquet") as fh:
        pf = pq.ParquetFile(fh)
        for i in range(pf.num_row_groups):
            t = pf.read_row_group(i, columns=["cell_name", "drug", "drugname_drugconc"])
            chunk = t.to_pandas()
            part = chunk.groupby(["cell_name", "drug", "drugname_drugconc"], sort=False).size()
            counts = part if counts is None else counts.add(part, fill_value=0)
            if verbose and i % 10 == 0:
                print(f"[tahoe] row group {i + 1}/{pf.num_row_groups}", flush=True)
    df = counts.astype("int64").rename("n_cells").reset_index()
    drug = df["drug"].astype(str).str.strip()
    vehicle = drug.isin(VEHICLES)
    df["intervention_kind"] = "chemical"
    df.loc[vehicle, "intervention_kind"] = "none"
    df["intervention_target"] = drug.where(~vehicle, "")
    df["intervention_detail"] = df["drugname_drugconc"].astype(str).map(_dose).where(~vehicle, "")
    df["identity"] = df["cell_name"].astype(str).map(normalise_cell_line)
    df["cell_line_id"] = df["cell_name"].map(cellosaurus).fillna("")
    df["source"] = "tahoe-100m"
    df["species"] = "human"
    df["identity_kind"] = "cell_line"
    df["assay"] = "Parse (Mosaic multiplexing)"
    df["licence"] = LICENCE
    return df
