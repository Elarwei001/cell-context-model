"""CZ CELLxGENE Census: observational single-cell contexts (metadata only).

Streams the ``obs`` table of each species (primary cells only, so cells re-used across collections
are counted once) and aggregates cell counts per (dataset, donor, tissue, cell type,
developmental stage, disease, assay).
"""

from __future__ import annotations

import re

import pandas as pd

SPECIES = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse",
    "macaca_mulatta": "rhesus macaque",
    "callithrix_jacchus": "common marmoset",
    "pan_troglodytes": "chimpanzee",
}
OBS_COLUMNS = [
    "dataset_id",
    "donor_id",
    "tissue",
    "tissue_general",
    "cell_type",
    "development_stage",
    "disease",
    "assay",
]
_SPATIAL = re.compile(r"visium|slide-seq|merfish|xenium|stereo|cosmx|seqfish|starmap", re.I)
LICENCE = "CC-BY-4.0 (per dataset; see CELLxGENE)"


def census_contexts(
    species: str = "homo_sapiens", census_version: str = "2025-11-08", verbose: bool = True
) -> pd.DataFrame:
    import cellxgene_census

    counts: pd.Series | None = None
    n = 0
    with cellxgene_census.open_soma(census_version=census_version) as census:
        obs = census["census_data"][species].obs
        reader = obs.read(column_names=OBS_COLUMNS, value_filter="is_primary_data == True")
        for table in reader:
            chunk = table.to_pandas()
            for col in OBS_COLUMNS:
                chunk[col] = chunk[col].astype(str)
            part = chunk.groupby(OBS_COLUMNS, sort=False).size()
            counts = part if counts is None else counts.add(part, fill_value=0)
            n += len(chunk)
            if verbose:
                print(f"[census:{species}] {n:,} cells, {len(counts):,} groups", flush=True)
    df = counts.astype("int64").rename("n_cells").reset_index()
    df["source"] = f"census:{census_version}:" + df["dataset_id"]
    df["modality"] = "single-cell RNA"
    df["species"] = SPECIES.get(species, species)
    df["identity_kind"] = "cell_type"
    df["identity"] = df["cell_type"]
    df["donor_id"] = df["dataset_id"] + ":" + df["donor_id"]
    df["intervention_kind"] = "none"
    df["has_spatial"] = df["assay"].str.contains(_SPATIAL)
    df["licence"] = LICENCE
    return df
