# What `gget` can already access (measured)

Measured on 2026-10-04 with [gget](https://github.com/scverse/gget) 0.30.8 and `cellxgene-census` 1.18.0 in an isolated Python 3.12 environment. Where gget wraps an upstream API, the corresponding public API was queried directly for counts.

## Summary

| Module | Source | What it returns | Diversity (measured) | Direction information | Value for this project |
|---|---|---|---|---|---|
| `cellxgene` | CZ CELLxGENE Census (stable = LTS 2025-11-08) | Single-cell expression and metadata, filterable by tissue, cell type, disease, donor, developmental stage | See below | None (observational); some datasets carry developmental stage or spatial coordinates | **Highest**: backbone for context diversity |
| `archs4` | ARCHS4 (uniformly processed GEO bulk RNA-seq) | Per-gene expression distribution across tissues and cell categories; gene correlations | Example: TP53 returns 72 tissue/cell categories | Summary statistics only; the underlying GEO samples include many treatments and time courses but need metadata mining | Medium: very broad contexts, but gget exposes summaries only |
| `bgee` | Bgee | Anatomical entities where a gene is expressed, with expression scores | Example: TP53 returns the top 100 anatomical entities | None | Low–medium: multi-species, multi-stage expression prior |
| `8cube` | 8cubeDB | snRNA-seq of 8 mouse strains × tissues × 4 individuals per sex; summary statistics (mean, variance, specificity) | Genetic-background diversity (8 strains) | Natural genetic variation (weak) | Low–medium: summaries only, no cell-level data |
| `cbio` | cBioPortal | Tumour and cell-line omics | 548 studies, about 414,000 samples; CCLE 2025: 1,981 cell lines | Mutation → expression (natural experiment, weak) | Medium: cell-line context diversity |
| `opentargets` | Open Targets | Gene–disease evidence, baseline expression, DepMap, interactions | — | Interactions partly directed | Low: prior |
| `enrichr` | Enrichr | Gene-set enrichment, including some pathway libraries | — | Pathway membership, no direction | Low: prior or evaluation |

Note: `gget cbio` search needs the optional `bravado` dependency (`gget setup cbio`); without it, searches return empty results.

## CELLxGENE Census (LTS 2025-11-08)

From `census_info/summary_cell_counts` (categories with at least one cell):

| Species | Cells (total / unique) | Cell types | Tissues (fine / general) | Diseases | Assays |
|---|---|---|---|---|---|
| Human | 162,025,130 / 99,633,637 | 903 | 423 / 71 | 261 | 39 |
| Mouse | 46,299,127 / 21,029,771 | 492 | 102 / 36 | 18 | 18 |
| Rhesus macaque | 7,010,229 / 2,929,014 | 54 | 29 / 2 | 1 | 2 |
| Common marmoset | 2,275,451 / 1,712,738 | 40 | 33 / 1 | 1 | 1 |
| Chimpanzee | 158,099 / 158,099 | 25 | 1 / 1 | 1 | 1 |

The Census contains 1,845 datasets from 313 collections.

## CELLxGENE Discover dataset metadata

The Discover curation API (`/curation/v1/datasets`, 2,237 datasets) covers more than the Census, including spatial datasets.

| | Human | Mouse |
|---|---|---|
| Datasets | 1,630 | 498 |
| Cells | about 210.8 million | about 69.1 million |
| Distinct developmental stages | 197 (41 embryonic or fetal) | 82 (16 Theiler stages) |
| Datasets spanning at least 3 developmental stages | 986 | 118 |
| Datasets with a spatial assay | 399 (Visium 319, Slide-seqV2 80, …) | 274 (mostly Slide-seqV2) |
| Donors (distinct donor IDs summed over datasets) | about 41,900 | — |

**Perturbations are almost absent.** Only 7 of 2,237 datasets carry perturbation metadata (`perturbation_types`): protein 5, chemical 1, diet 1. None declares a `genetic_perturbation_strategy`. Observational atlases therefore supply breadth of contexts but no direction; direction must come from other sources (see [open-data-survey.md](open-data-survey.md)).

## Reproducing

```python
import cellxgene_census

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    summary = census["census_info"]["summary_cell_counts"].read().concat().to_pandas()
    datasets = census["census_info"]["datasets"].read().concat().to_pandas()
```

Discover metadata: `GET https://api.cellxgene.cziscience.com/curation/v1/datasets` (about 10 MB of JSON).
