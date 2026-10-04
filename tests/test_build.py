import pandas as pd

from ccm.catalog.build import build_contexts, build_identities, build_interventions
from ccm.catalog.schema import normalise_cell_line


def _rows(**common):
    def row(**kw):
        r = dict(common)
        r.update(kw)
        return r

    return row


def test_normalise_cell_line():
    assert normalise_cell_line("Hep G2") == "HEPG2"
    assert normalise_cell_line("HCT-116") == "HCT116"
    assert normalise_cell_line("hTERT RPE-1") == "RPE1"


def test_contexts_merge_across_sources_and_count_cells():
    line = _rows(species="human", identity_kind="cell_line", intervention_kind="genetic")
    a = pd.DataFrame(
        [
            line(source="s1", identity="K562", intervention_target="GATA1", n_cells=100),
            line(source="s1", identity="K562", intervention_kind="none", n_cells=1000),
        ]
    )
    b = pd.DataFrame([line(source="s2", identity="K562", intervention_target="GATA1", n_cells=50)])
    ctx = build_contexts([a, b])
    gata1 = ctx[ctx["intervention_target"] == "GATA1"].iloc[0]
    assert gata1["n_cells"] == 150
    assert gata1["n_sources"] == 2
    assert gata1["sources"] == "s1,s2"
    assert gata1["direction_kinds"] == "genetic"
    assert ctx["context_id"].is_unique


def test_interventions_count_distinct_cell_lines_and_cell_types():
    rows = pd.DataFrame(
        [
            dict(source="p1", species="human", identity_kind="cell_line", identity=line,
                 intervention_kind="genetic", intervention_target="MYC", n_cells=10)
            for line in ("K562", "RPE1", "JURKAT")
        ]
        + [
            dict(source="p2", species="human", identity_kind="cell_type", identity="CD4 T cell",
                 donor_id=d, intervention_kind="genetic", intervention_target="MYC", n_cells=5)
            for d in ("d1", "d2")
        ]
        + [
            dict(source="p1", species="human", identity_kind="cell_line", identity="K562",
                 intervention_kind="genetic", intervention_target="GATA1", n_cells=7)
        ]
    )
    iv = build_interventions(build_contexts([rows])).set_index("intervention_target")
    assert iv.loc["MYC", "n_cell_lines"] == 3
    assert iv.loc["MYC", "n_cell_types"] == 1
    assert iv.loc["MYC", "n_identities"] == 4
    assert iv.loc["MYC", "n_donors"] == 2
    assert iv.loc["GATA1", "n_identities"] == 1
    assert iv.index[0] == "MYC"  # sorted by identity coverage


def test_identities_split_observational_and_perturbed():
    rows = pd.DataFrame(
        [
            dict(source="c", species="human", identity_kind="cell_type", identity="hepatocyte",
                 tissue="liver", development_stage=stage, n_cells=100)
            for stage in ("fetal 10w", "fetal 12w", "adult")
        ]
        + [
            dict(source="p", species="human", identity_kind="cell_type", identity="hepatocyte",
                 intervention_kind="chemical", intervention_target="drugX", n_cells=30),
            dict(source="p", species="human", identity_kind="cell_type", identity="hepatocyte",
                 intervention_kind="genetic", intervention_target="HNF4A", n_cells=20),
        ]
    )
    ids = build_identities(build_contexts([rows])).set_index("identity")
    h = ids.loc["hepatocyte"]
    assert h["n_cells_observational"] == 300
    assert h["n_cells_perturbed"] == 50
    assert h["n_genetic"] == 1 and h["n_chemical"] == 1 and h["n_interventions"] == 2
    assert h["n_development_stages"] == 3
    assert bool(h["has_temporal"])
