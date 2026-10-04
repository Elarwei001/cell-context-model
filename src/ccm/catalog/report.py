"""Render catalog/summary.md and catalog/gaps.md from the three catalog tables."""

from __future__ import annotations

import numpy as np
import pandas as pd

THRESHOLDS = (1, 2, 5, 10, 20, 50, 100)


def _md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(_fmt(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def _fmt(v) -> str:
    if isinstance(v, (int,)) or (hasattr(v, "dtype") and "int" in str(getattr(v, "dtype", ""))):
        return f"{int(v):,}"
    if isinstance(v, float):
        return f"{v:,.3g}"
    return str(v)


def identity_counts(
    contexts: pd.DataFrame, by: str | None, per_source: dict[str, pd.DataFrame] | None = None
) -> pd.DataFrame:
    """Distinct cell types, cell lines, tissues, donors, contexts and interventions.

    ``by`` is None (overall), "species", or "source". For "source", ``per_source`` maps each
    source name to the contexts built from that source alone, so shared contexts are counted
    in every source that contributes them.
    """
    if by == "source":
        groups = list((per_source or {}).items())
    elif by in ("species", "modality"):
        groups = list(contexts.groupby(by))
    else:
        groups = [("all", contexts)]
    rows = []
    for name, g in groups:
        pert = g[g["intervention_kind"] != "none"]
        rows.append(
            {
                by or "scope": name,
                "cell types": g.loc[g["identity_kind"] == "cell_type", "identity"].nunique(),
                "cell lines": g.loc[g["identity_kind"] == "cell_line", "identity"].nunique(),
                "tissues": g.loc[g["tissue"] != "", "tissue"].nunique(),
                "donors": g.loc[g["donor_id"] != "", "donor_id"].nunique(),
                "contexts": len(g),
                "interventions": pert[["intervention_kind", "intervention_target"]]
                .drop_duplicates()
                .shape[0],
                "identities with interventions": pert[["identity_kind", "identity"]]
                .drop_duplicates()
                .shape[0],
                "cells": int(g["n_cells"].sum()),
                "bulk samples": int(g["n_samples"].sum()),
            }
        )
    return pd.DataFrame(rows)


def coverage_distribution(interventions: pd.DataFrame) -> pd.DataFrame:
    """Number of interventions measured in at least k distinct identities, per kind."""
    rows = []
    if interventions.empty:
        return pd.DataFrame()
    keyed = interventions.assign(
        _k=interventions["intervention_kind"]
        + np.where(interventions["intervention_effect"] != "",
                   " (" + interventions["intervention_effect"] + ")", "")
    )
    for kind, g in keyed.groupby("_k"):
        row = {"intervention kind": kind, "interventions": len(g)}
        for k in THRESHOLDS:
            row[f">= {k} identities"] = int((g["n_identities"] >= k).sum())
        row["max identities"] = int(g["n_identities"].max())
        rows.append(row)
    return pd.DataFrame(rows)


def concentration(contexts: pd.DataFrame) -> pd.DataFrame:
    """How concentrated cells are over contexts (single-cell contexts only)."""
    sc = contexts.loc[contexts["n_cells"] > 0, "n_cells"].sort_values(ascending=False)
    total = sc.sum()
    rows = []
    for frac in (0.001, 0.01, 0.1):
        k = max(1, int(round(len(sc) * frac)))
        rows.append({"top share of contexts": f"{frac:.1%}", "contexts": k,
                     "share of cells": f"{sc.iloc[:k].sum() / total:.1%}"})
    rows.append({"top share of contexts": "median cells per context", "contexts": len(sc),
                 "share of cells": f"{int(sc.median()):,}"})
    return pd.DataFrame(rows)


def direction_coverage(identities: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for kind, g in identities.groupby("identity_kind"):
        rows.append(
            {
                "identity kind": kind,
                "identities": len(g),
                "with genetic": int((g["n_genetic"] > 0).sum()),
                "with chemical": int((g["n_chemical"] > 0).sum()),
                "with ligand": int((g["n_ligand"] > 0).sum()),
                "with >=2 stages or time points": int(g["has_temporal"].sum()),
                "with imaging": int(g["has_imaging"].sum()),
                "with spatial": int(g["has_spatial"].sum()),
                "observational only": int((g["n_interventions"] == 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def render_summary(
    contexts, interventions, identities, sources_note: str, per_source=None
) -> str:
    from .build import build_interventions

    by_mod = {
        m: build_interventions(contexts[contexts["modality"] == m])
        for m in sorted(contexts["modality"].unique())
    }
    sc_iv = by_mod.get("single-cell RNA", build_interventions(contexts.iloc[:0]))
    cols = ["intervention_target", "n_identities", "n_cell_lines", "n_cell_types", "n_cells",
            "n_samples", "identities"]
    def _loss(iv):
        return iv[(iv["intervention_kind"] == "genetic") & (iv["intervention_effect"] == "loss")]

    top_sc = _loss(sc_iv).head(15)[cols]
    top_any = _loss(interventions).head(10)[cols[:-1] + ["modalities"]]
    parts = [
        "# Context catalog: summary",
        "",
        "Generated by `scripts/build_catalog.py`. " + sources_note,
        "",
        "## Identity counts",
        "",
        "### Overall",
        "",
        _md_table(identity_counts(contexts, None)),
        "",
        "### By species",
        "",
        _md_table(identity_counts(contexts, "species")),
        "",
        "### By modality",
        "",
        _md_table(identity_counts(contexts, "modality")),
        "",
        "### By source",
        "",
        _md_table(identity_counts(contexts, "source", per_source)),
        "",
        "## In how many distinct cell lines / cell types was each intervention measured?",
        "",
        "Each cell counts the interventions measured in at least *k* distinct identities "
        "(cell lines + cell types), separately for each measurement modality (bulk L1000 "
        "signatures cover many more cell lines but only about 978 genes; imaging counts wells), "
        "then for all modalities combined.",
        "",
        *[
            line
            for m, iv in by_mod.items()
            if not iv.empty
            for line in (f"### {m}", "", _md_table(coverage_distribution(iv)), "")
        ],
        "### All data combined",
        "",
        _md_table(coverage_distribution(interventions)),
        "",
        "### Loss-of-function genetic perturbations with the widest single-cell identity coverage",
        "",
        _md_table(top_sc),
        "",
        "### Loss-of-function genetic perturbations with the widest coverage (any modality)",
        "",
        _md_table(top_any),
        "",
        "Note: for cell types, \"stages or time points\" mostly means observational cells "
        "sampled at two or more developmental stages, not time-resolved interventions.",
        "",
        "## Direction-bearing data per identity",
        "",
        _md_table(direction_coverage(identities)),
        "",
        "## Concentration of cells over contexts",
        "",
        _md_table(concentration(contexts)),
        "",
    ]
    return "\n".join(parts)


def render_gaps(interventions, identities) -> str:
    single = interventions[interventions["n_identities"] == 1]
    by_kind = single.groupby("intervention_kind").size().rename("measured in one identity only")
    obs_only = identities[identities["n_interventions"] == 0]
    big_obs_only = obs_only.sort_values("n_cells_observational", ascending=False).head(25)[
        ["identity_kind", "identity", "species", "n_cells_observational", "n_tissues",
         "n_development_stages"]
    ]
    parts = [
        "# Context catalog: gaps",
        "",
        "## Interventions measured in a single identity",
        "",
        _md_table(by_kind.reset_index()),
        "",
        f"## Identities with observational data only: {len(obs_only):,} of {len(identities):,}",
        "",
        "The 25 largest by number of observational cells:",
        "",
        _md_table(big_obs_only),
        "",
    ]
    return "\n".join(parts)
