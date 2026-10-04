# Context catalog: design

Status: design draft, 2026-10-04.

## 1. What it is

The context catalog is a single table that regroups cells from many open datasets by **the context they are in**, with one row per context.

A *context* is the combination of conditions that determines a cell's state:

| Axis | Examples |
|---|---|
| Species | human, mouse, macaque |
| Biological identity | tissue and cell type (hepatocyte in liver, CD4 T cell in blood), or a cell line (K562, A549) |
| State | developmental stage, disease, donor, sex, age, time after stimulation |
| Intervention | none; a genetic perturbation (CRISPRi/CRISPRko of gene X); a chemical (drug Y at dose d); a ligand or cytokine (IFN-γ) |

Example rows (illustrative values):

| species | identity | state | intervention | n_cells | sources | direction information |
|---|---|---|---|---|---|---|
| human | CD4 T cell | resting, donor 1 | CRISPRi IL2RA | 180 | primary CD4 T Perturb-seq | genetic intervention |
| human | CD4 T cell | stimulated 8 h, donor 1 | none | 52,000 | primary CD4 T Perturb-seq (controls) | time course |
| human | CD4 T cell | healthy adult | none | 1,200,000 | 37 CELLxGENE immune datasets | none |
| human | PBMC monocyte | donor 3 | IFN-γ | 9,000 | Parse 10M PBMC | ligand intervention |
| human | hepatocyte | fetal, 12 weeks post-conception | none | 30,000 | CELLxGENE fetal liver atlas | developmental stage |
| human | A549 | — | drug X | 3,000 | Tahoe-100M | chemical intervention |

The catalog is built from **metadata only**. No expression matrices are downloaded, so it is cheap to build and to refresh.

## 2. Why build it first

1. **Measure diversity honestly.** "100 million cells" says little about diversity if most cells sit in a few contexts. Counting contexts, and how cells are distributed over them, is the right measure for a project that values diversity over cell count.
2. **Drive pretraining sampling.** Cap the number of cells drawn per context (for example 2,000) so that the model sees as many distinct contexts as possible and is not dominated by a few very large datasets.
3. **Count identities, not just cells.** For every intervention, report in how many distinct cell lines and cell types it has been measured; for every cell line or cell type, report how many interventions and which kinds of direction-bearing data exist (section 5). Learning how an effect depends on context requires the same intervention across many identities.
4. **Find the gaps.** Which identities have only observational data and no intervention? Which interventions exist in only one cell line? Which tissues lack time courses? The gaps decide what data to look for next.
5. **Define evaluation.** Held-out splits can be expressed directly as sets of contexts, such as a whole cell type, a donor or a cell line. That makes "generalisation to unseen contexts" testable.

## 3. Schema (draft)

One row per context key; the context key is the tuple of normalised axis values.

| Column | Type | Notes |
|---|---|---|
| `context_id` | string | Stable hash of the context key |
| `species` | string | NCBI taxon label |
| `identity_kind` | enum | `cell_type`, `cell_line`, `organoid`, `mixed` |
| `tissue`, `tissue_general` | string | UBERON labels where available |
| `cell_type` | string | Cell Ontology (CL) label where available |
| `cell_line` | string | Cellosaurus ID where available |
| `development_stage` | string | HsapDv/MmusDv label where available |
| `disease` | string | MONDO label; `normal` for healthy |
| `donor_id` | string | Source-scoped donor identifier (not linked across sources) |
| `time_point` | string | Time after stimulation or intervention, if any |
| `intervention_kind` | enum | `none`, `genetic`, `chemical`, `ligand`, `other` |
| `intervention_target` | string | Gene symbol (genetic), compound (chemical) or ligand name |
| `intervention_detail` | string | Modality (CRISPRi/CRISPRko/overexpression), dose, guide, etc. |
| `n_cells` | int | Cells in this context, summed over sources (primary data only where the source marks duplicates) |
| `n_sources` | int | Number of contributing datasets |
| `sources` | list | Dataset identifiers |
| `assays` | list | Assay labels (10x 3' v3, Smart-seq2, Visium, ...) |
| `has_spatial` | bool | Spatial coordinates available |
| `has_kinetics` | bool | Spliced/unspliced counts or metabolic labelling available |
| `direction_kinds` | list | Subset of {`genetic`, `chemical`, `ligand`, `temporal`, `lineage`, `spatial`, `natural_variation`} |
| `licence` | string | Most restrictive licence among contributing sources |

Ontology normalisation reuses the source annotations: CELLxGENE already ships CL, UBERON, MONDO and HsapDv/MmusDv terms. For sources without ontology terms, a mapping table is kept in the repository and reviewed by hand.

## 4. First sources (phase 1)

| Source | Access | What it contributes |
|---|---|---|
| CELLxGENE Census LTS 2025-11-08 | `cellxgene_census` / `gget cellxgene` (obs metadata only) | Observational breadth: 903 human and 492 mouse cell types; developmental stages; spatial datasets via the Discover API |
| Arc scBaseCount | Sample-level metadata (CC0) | Additional observational breadth across 27 species; kinetics-ready counts |
| Primary CD4 T genome-wide Perturb-seq | Per-cell metadata | Genetic interventions × donors × stimulation time points |
| Replogle 2022, Nadig 2025, X-Atlas/Orion | Per-cell metadata | Genetic interventions across six cell lines |
| Parse 10M PBMC | Per-cell metadata (registration required) | Ligand interventions × donors |
| Tahoe-100M | Per-cell metadata | Chemical interventions × 50 cell lines |
| LINCS L1000 | Signature metadata | Genetic and chemical interventions × about 250 cell lines (bulk) |

## 5. Identity coverage: how many cell lines and cell types

The context table alone does not show the number that matters most for learning context dependence: **across how many distinct cell lines or cell types** a given kind of data, or a given intervention, has been measured. A knockdown measured in 1 cell line teaches nothing about how its effect changes with context; the same knockdown in 50 cell lines does. The catalog therefore ships two derived views keyed on identity counts.

### 5.1 Intervention coverage (`catalog/interventions.parquet`)

One row per intervention (for example "CRISPRi IL2RA", "IFN-γ", "drug X").

| Column | Notes |
|---|---|
| `intervention_kind`, `intervention_target` | As in the context table |
| `n_cell_lines` | Distinct cell lines in which this intervention was measured |
| `n_cell_types` | Distinct primary cell types (CL terms) in which it was measured |
| `n_identities` | `n_cell_lines + n_cell_types` (identities are counted once) |
| `n_donors`, `n_states` | Distinct donors and states (time points, stimulation) |
| `n_cells` | Total perturbed cells |
| `identities` | List of the cell lines and cell types |
| `sources` | Contributing datasets |

The headline statistic is the distribution of `n_identities` over interventions. For example: how many genes have been knocked down in at least 5 distinct cell lines or cell types? In at least 20?

### 5.2 Identity coverage (`catalog/identities.parquet`)

One row per identity (cell line or cell type).

| Column | Notes |
|---|---|
| `identity_kind`, `identity` | Cell line (Cellosaurus) or cell type (CL) |
| `n_contexts` | Contexts for this identity |
| `n_states`, `n_donors` | Distinct states and donors |
| `n_genetic`, `n_chemical`, `n_ligand` | Distinct interventions of each kind measured in this identity |
| `has_temporal`, `has_spatial`, `has_kinetics` | Direction-bearing data available for this identity |
| `n_cells_observational`, `n_cells_perturbed` | Cell counts |

The headline statistics:
- how many distinct cell lines and cell types exist in total;
- how many of them have any intervention data;
- how many have each kind of direction-bearing data.

### 5.3 Catalog-level counts

Every catalog build reports, per species and per source and overall:
- distinct cell types;
- distinct cell lines;
- distinct tissues;
- distinct donors;
- distinct contexts;
- distinct interventions;
- the number of identities with at least one intervention.

## 6. Outputs

- `catalog/contexts.parquet`: the context table (section 3).
- `catalog/interventions.parquet`: intervention coverage (section 5.1).
- `catalog/identities.parquet`: identity coverage (section 5.2).
- `catalog/summary.md`: generated summary, containing:
  - the catalog-level identity counts (section 5.3);
  - the distribution of `n_identities` per intervention;
  - contexts per axis;
  - the cell-count distribution over contexts (for example the share of cells held by the top 1% of contexts);
  - coverage by direction kind.
- `catalog/gaps.md`: generated list of identity × direction-kind combinations with no data.

## 7. Open questions

- **Granularity.** How fine should the context key be? Should donor be part of the key for observational data, or only for perturbation data? A key that is too fine fragments the catalog; a key that is too coarse hides diversity. The initial proposal is to include donor and keep a coarser `context_group` column for sampling.
- **Cross-source identity.** The same cell line or cell type can be labelled differently across sources. CL and Cellosaurus IDs resolve most cases; the remainder go into a reviewed mapping table.
- **Duplicates.** Cells reused across collections should be counted once. CELLxGENE marks primary data; other sources need source-level deduplication rules.
