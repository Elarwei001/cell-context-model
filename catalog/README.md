# Context catalog

Generated tables describing which cellular contexts open data covers. Design: [../docs/design/context-catalog.md](../docs/design/context-catalog.md).

| File | Committed | Content |
|---|---|---|
| `summary.md` | yes | Identity counts (cell types, cell lines, tissues, donors, contexts, interventions), overall, by species and by source. How many interventions were measured in ≥ k distinct identities. Direction-bearing data per identity. Concentration of cells over contexts |
| `gaps.md` | yes | Interventions measured in a single identity; identities with observational data only |
| `contexts.parquet` | no (large) | One row per context |
| `interventions.parquet` | no | One row per intervention, with `n_cell_lines`, `n_cell_types`, `n_identities`, `n_donors`, `n_states` |
| `identities.parquet` | no | One row per cell line or cell type, with intervention and direction coverage |

## Rebuilding

```bash
uv venv -p 3.12 .venv
VIRTUAL_ENV=.venv uv pip install -e ".[census,h5ad,hf]"
.venv/bin/python scripts/build_catalog.py --data-dir data --sources census,h5ad,xatlas,tahoe,l1000
```

- `census`, `xatlas` and `tahoe` are read remotely (CELLxGENE Census on S3; Hugging Face with column projection). Only metadata columns are transferred.
- `l1000` downloads two text files from the CMap LINCS 2020 release (about 465 MB).
- `h5ad` expects the single-cell CRISPR screens as local h5ad files in `data/h5ad/` (file names in `scripts/build_catalog.py`):
  - Replogle 2022 K562 genome-wide and RPE1, via scPerturb (Zenodo 7041849);
  - Replogle 2022 K562 essential (figshare);
  - Nadig 2025 Jurkat and HepG2 (GEO GSE264667).

  Only `obs` is read.

## Caveats

- **Cell lines** are matched across sources by a normalised name (upper case, alphanumerics only, a few aliases) and, where available, a Cellosaurus accession. Rare naming differences can still split one line into two identities.
- **Cell types** are CELLxGENE Cell Ontology labels. Primary-cell perturbation datasets with cell-type identities are not yet included, so every intervention row currently comes from cell lines.
- **Identities are pooled across species** in `identities.parquet` and in `gaps.md`; for example, `oligodendrocyte` from all five species is one row. Per-species counts are in the "By species" table of `summary.md`.
- **Donors** are scoped to their source dataset and never linked across datasets.
- **L1000** contributes bulk signatures (`n_samples`), not single cells. Its shRNA seed-matched controls (`trt_sh.css`) and antibody treatments (`trt_aby`) are counted as `other`.
