# cell-context-model

Learning how genes and cells influence each other across diverse cellular contexts, from open single-cell and perturbation data.

## Goal

We want to train a transformer-type model that captures two kinds of biological knowledge:

1. **Within-cell regulation.** When an upstream gene changes, which downstream genes respond, and in which direction?
2. **Between-cell influence.** How do signals from one cell (for example a ligand or a cytokine) change the expression of another cell?

A model that has learned these relations should be able to predict the effect of a perturbation in a cellular context it has never seen perturbed. For example, it should predict the response to knocking down a gene in a new cell type, given only that cell type's unperturbed state.

## Guiding principles

- **Diversity of contexts over number of cells.** A context is the set of conditions that defines a cell's state: species, tissue, cell type, developmental stage, disease, donor, stimulation, and intervention. One hundred million cells concentrated in a few contexts teach less than a few million cells spread across thousands of contexts. Recent evidence supports this; see [docs/literature/2026-literature-review.md](docs/literature/2026-literature-review.md).
- **Direction needs more than observation.** Observational atlases show which genes co-vary, not which gene drives which. Direction comes from:
  - **interventions**: genetic, chemical, and ligand perturbations;
  - **temporal order**: development, differentiation, stimulation time courses, and RNA kinetics;
  - **spatial neighbourhoods**: cell–cell communication.
- **Evaluate across contexts from day one.** Held-out evaluation splits are defined by unseen contexts (cell lines, donors, tissues, conditions), never by random cells.
- **Open data, tracked licences.** Each data source is recorded with its licence. Non-commercial sources are kept separate from permissive ones.

## Status

Early stage: data inventory and design.

| Area | Document |
|---|---|
| What the `gget` library can already access, measured on 2026-10-04 | [docs/data/gget-coverage.md](docs/data/gget-coverage.md) |
| Survey of about 45 open data sources beyond gget, ranked by context diversity and direction information | [docs/data/open-data-survey.md](docs/data/open-data-survey.md) |
| Design of the **context catalog**, the first deliverable | [docs/design/context-catalog.md](docs/design/context-catalog.md) |
| 2026 literature review: cross-context perturbation prediction, in-context learning over cell sets, flow matching, knowledge priors, benchmarks | [docs/literature/2026-literature-review.md](docs/literature/2026-literature-review.md) |
| Deep dive into Stack (Arc Institute), the closest existing in-context single-cell model | [docs/literature/stack-deep-dive.md](docs/literature/stack-deep-dive.md) |
| Close reading of CellMSA and PT-RAG | [docs/literature/cellmsa-and-pt-rag.md](docs/literature/cellmsa-and-pt-rag.md) |

## Roadmap

1. **Context catalog.** Merge metadata from CELLxGENE Census, scBaseCount, and the major perturbation datasets into a single table, with one row per context. Record:
   - the number of cells in each context;
   - which kinds of direction information each context has;
   - the number of distinct cell lines and cell types: overall, per source, and for each intervention (in how many identities it has been measured).

   This shows how diverse the open data really is, drives pretraining sampling, and exposes the gaps.
2. **Pilot domain.** Pick one domain that has observational data, genetic perturbations, ligand perturbations and time points at a manageable scale; primary human CD4 T cells are the leading candidate. Use it to test whether a model learns upstream/downstream relations that transfer across donors and stimulation states.
3. **Model.** Design and train a model that combines a diverse observational corpus, direction-bearing data, and perturbation data, evaluated on held-out contexts.
4. **Tooling.** Contribute missing data-access modules (for example directed signalling networks and a perturbation-dataset catalog) to [gget](https://github.com/scverse/gget).

## Licence

Code and documentation in this repository are released under the [Apache License 2.0](LICENSE). The licences of the third-party datasets and models discussed here are their own. Several are non-commercial; see the per-source notes in [docs/data/open-data-survey.md](docs/data/open-data-survey.md).
