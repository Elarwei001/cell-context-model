"""Single-cell CRISPR screens distributed as h5ad (Replogle 2022 via scPerturb, Nadig 2025 / GEO).

Only ``obs`` is read (backed mode), so multi-GB files are cheap to scan. One context row is
emitted per target gene, plus one for the non-targeting controls (``intervention_kind="none"``).
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from ..schema import CELLOSAURUS, normalise_cell_line

CONTROL_LABELS = {"non-targeting", "control", "non_targeting", "ntc"}


@dataclass(frozen=True)
class H5adScreen:
    source: str  # e.g. "replogle2022:k562_gwps"
    path: str
    cell_line: str
    modality: str = "CRISPRi"
    target_column: str = "gene"
    licence: str = ""
    assay: str = "10x 3'"


def screen_contexts(screen: H5adScreen) -> pd.DataFrame:
    import anndata as ad

    adata = ad.read_h5ad(screen.path, backed="r")
    targets = adata.obs[screen.target_column].astype(str)
    adata.file.close()
    counts = targets.value_counts()
    line = normalise_cell_line(screen.cell_line)
    df = pd.DataFrame({"intervention_target": counts.index, "n_cells": counts.to_numpy()})
    is_ctrl = df["intervention_target"].str.lower().isin(CONTROL_LABELS)
    df["intervention_kind"] = "genetic"
    df.loc[is_ctrl, "intervention_kind"] = "none"
    df.loc[is_ctrl, "intervention_target"] = ""
    df = df.groupby(["intervention_kind", "intervention_target"], as_index=False)["n_cells"].sum()
    df["intervention_detail"] = (df["intervention_kind"] == "genetic").map(
        {True: screen.modality, False: ""}
    )
    df["source"] = screen.source
    df["modality"] = "single-cell RNA"
    df["species"] = "human"
    df["identity_kind"] = "cell_line"
    df["identity"] = line
    df["cell_line_id"] = CELLOSAURUS.get(line, "")
    df["assay"] = screen.assay
    df["licence"] = screen.licence
    return df
