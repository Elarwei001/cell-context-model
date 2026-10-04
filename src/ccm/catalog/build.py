"""Merge source rows into the three catalog tables and render the summary and gap reports.

Tables (see docs/design/context-catalog.md):
- contexts: one row per context key, summed over sources;
- interventions: one row per (intervention_kind, intervention_target), with the number of
  distinct cell lines / cell types / donors / states in which it was measured;
- identities: one row per (identity_kind, identity), with intervention and direction coverage.
"""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

from .schema import KEY_COLUMNS, conform

STATE_COLUMNS = ["tissue", "development_stage", "disease", "time_point"]


def _context_id(row: pd.Series) -> str:
    raw = "\x1f".join(str(row[c]) for c in KEY_COLUMNS)
    return hashlib.sha1(raw.encode()).hexdigest()[:16]


def _direction_kinds(row: pd.Series) -> str:
    kinds = []
    if row["intervention_kind"] != "none":
        kinds.append(row["intervention_kind"])
    if row["time_point"]:
        kinds.append("temporal")
    if row["has_spatial"]:
        kinds.append("spatial")
    if row["has_kinetics"]:
        kinds.append("kinetics")
    return ",".join(kinds)


def _join_unique(values) -> str:
    return ",".join(sorted({v for v in values if v}))


def build_contexts(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Merge per-source rows into one row per context key."""
    rows = pd.concat([conform(f) for f in frames], ignore_index=True)
    agg = (
        rows.groupby(KEY_COLUMNS, sort=False, dropna=False)
        .agg(
            cell_line_id=("cell_line_id", _join_unique),
            tissue_general=("tissue_general", _join_unique),
            cell_type=("cell_type", _join_unique),
            n_cells=("n_cells", "sum"),
            n_samples=("n_samples", "sum"),
            n_sources=("source", "nunique"),
            sources=("source", _join_unique),
            assays=("assay", _join_unique),
            has_spatial=("has_spatial", "any"),
            has_kinetics=("has_kinetics", "any"),
            licence=("licence", _join_unique),
        )
        .reset_index()
    )
    agg["direction_kinds"] = agg.apply(_direction_kinds, axis=1)
    agg.insert(0, "context_id", agg.apply(_context_id, axis=1))
    return agg


def _state_key(df: pd.DataFrame) -> pd.Series:
    return df[STATE_COLUMNS].astype(str).agg("|".join, axis=1)


def build_interventions(contexts: pd.DataFrame) -> pd.DataFrame:
    """One row per intervention, counting the distinct identities it was measured in."""
    pert = contexts[contexts["intervention_kind"] != "none"].copy()
    if pert.empty:
        return pd.DataFrame()
    pert["_state"] = _state_key(pert)
    pert["_line"] = np.where(pert["identity_kind"] == "cell_line", pert["identity"], "")
    pert["_type"] = np.where(pert["identity_kind"] == "cell_type", pert["identity"], "")
    pert["_donor"] = np.where(pert["donor_id"] != "", pert["sources"] + ":" + pert["donor_id"], "")
    g = pert.groupby(["intervention_kind", "intervention_target"], sort=False)
    out = g.agg(
        n_cell_lines=("_line", lambda s: s[s != ""].nunique()),
        n_cell_types=("_type", lambda s: s[s != ""].nunique()),
        n_donors=("_donor", lambda s: s[s != ""].nunique()),
        n_states=("_state", "nunique"),
        n_contexts=("context_id", "nunique"),
        n_cells=("n_cells", "sum"),
        n_samples=("n_samples", "sum"),
        species=("species", _join_unique),
        identities=("identity", _join_unique),
        sources=("sources", lambda s: _join_unique(",".join(s).split(","))),
    ).reset_index()
    out.insert(4, "n_identities", out["n_cell_lines"] + out["n_cell_types"])
    return out.sort_values(["n_identities", "n_cells"], ascending=False, ignore_index=True)


def build_identities(contexts: pd.DataFrame) -> pd.DataFrame:
    """One row per identity (cell line or cell type), with intervention and direction coverage."""
    c = contexts.copy()
    c["_state"] = _state_key(c)
    c["_donor"] = np.where(c["donor_id"] != "", c["sources"] + ":" + c["donor_id"], "")
    rows = []
    for (kind, ident), g in c.groupby(["identity_kind", "identity"], sort=False):
        pert = g[g["intervention_kind"] != "none"]
        obs = g[g["intervention_kind"] == "none"]
        rows.append(
            {
                "identity_kind": kind,
                "identity": ident,
                "species": _join_unique(g["species"]),
                "cell_line_id": _join_unique(g["cell_line_id"]),
                "n_contexts": len(g),
                "n_tissues": g.loc[g["tissue"] != "", "tissue"].nunique(),
                "n_development_stages": g.loc[
                    g["development_stage"] != "", "development_stage"
                ].nunique(),
                "n_states": g["_state"].nunique(),
                "n_donors": g.loc[g["_donor"] != "", "_donor"].nunique(),
                "n_genetic": pert.loc[pert["intervention_kind"] == "genetic", "intervention_target"]
                .nunique(),
                "n_chemical": pert.loc[
                    pert["intervention_kind"] == "chemical", "intervention_target"
                ].nunique(),
                "n_ligand": pert.loc[pert["intervention_kind"] == "ligand", "intervention_target"]
                .nunique(),
                "has_temporal": bool(
                    (g["time_point"] != "").any()
                    or g.loc[g["development_stage"] != "", "development_stage"].nunique() >= 2
                ),
                "has_spatial": bool(g["has_spatial"].any()),
                "has_kinetics": bool(g["has_kinetics"].any()),
                "n_cells_observational": int(obs["n_cells"].sum()),
                "n_cells_perturbed": int(pert["n_cells"].sum()),
                "n_samples": int(g["n_samples"].sum()),
                "sources": _join_unique(",".join(g["sources"]).split(",")),
            }
        )
    out = pd.DataFrame(rows)
    out["n_interventions"] = out["n_genetic"] + out["n_chemical"] + out["n_ligand"]
    return out.sort_values(
        ["n_interventions", "n_cells_observational"], ascending=False, ignore_index=True
    )
