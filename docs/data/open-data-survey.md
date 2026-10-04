# Open data for learning gene and cell upstream/downstream relations across diverse contexts

**Date:** 2026-10-04

**Scope.** This survey inventories open data for pretraining a transformer-type model that learns (a) within-cell regulatory direction (which genes act upstream or downstream of which) and (b) between-cell influence (how one cell's signals change another cell's state). The guiding principle is **diversity of cellular contexts/states over number of cells**.

**Verification level.** Compiled from web search, reading source pages, and live queries against the Hugging Face API and GitHub API. Only what was actually read is reported; any number not backed by a source read during this survey is marked **(unverified)**. Sizes used in the ranking sections are often order-of-magnitude estimates.

**gget module baseline** (`scverse/gget` main branch, checked 2026-10): 8cube, alphafold, archs4, bgee, blast, blat, cbio, cellxgene, cosmic, diamond, elm, enrichr, g2p, gpt, info, muscle, mutate, opentargets, pdb, ref, search, seq, setup, virus.

---

## 0. Current gget coverage (parts relevant to this project, confirmed from source code)

| gget module | What it actually returns | Relevance to this project |
|---|---|---|
| `cellxgene` | Pulls Census AnnData/metadata via `cellxgene_census.open_soma(census_version=...)`; default `"stable"`; a source comment notes that non-human primates require `census_version='2025-11-08'` (LTS) or newer | The only entry point covering large-scale observational single-cell data. Note gget issue #218 "Replace cellxgene-census" (upstream package maintenance / Python ≥3.14 compatibility; still OPEN) |
| `archs4` | ARCHS4 gene correlation (`which='correlation'`) and tissue expression (`'tissue'`) | Returns summary statistics only; does not download expression matrices |
| `bgee` | Orthologs/expression (Bgee API) | Cross-species expression background; no direction information |
| `8cube` | ψ/ζ tissue specificity and expression from 8cubeDB (`eightcubedb.onrender.com`) (mouse; e.g. `Akr1c21`) | Multi-strain × tissue background; no direction |
| `opentargets` | resource ∈ diseases/drugs/tractability/pharmacogenetics/expression/**depmap**/**interactions** (with a `sourceDatabase` field) | DepMap provides only gene-level gene-effect summaries; interactions are PPIs (the source-database field distinguishes origins; the exact set of sources is unverified) |
| `enrichr` | Enrichment analysis; shortcut libraries KEGG_2021_Human / ChEA_2016 / GO_BP_2021 / GWAS_Catalog_2019 / PanglaoDB_2021 / KEA_2015, or any Enrichr library name | Enrichment only; no edge tables returned. Whether LINCS/Reactome and similar library names are all available in Enrichr was not checked one by one (unverified) |
| `cbio`, `cosmic` | Tumor mutations / clinical data | Not core to this project |

Related gget PRs that were proposed recently but **closed without merging**: `reactome` (#236), `encode` (#231), `rummagene/rummageo` (#230), `alliance` (#232), MSigDB gene sets (#241). Issue `depmap` (#121) is closed (partly superseded by the opentargets depmap resource). In short: **gget currently has no dedicated module for perturbation data, spatial data, directed signaling networks, or eQTLs.**

---

## A. Large-scale observational single-cell corpora (the main source of context breadth)

| # | Name / URL | Organism | Modality | Scale | Diversity (source) | Direction information | Access / size | License | gget |
|---|---|---|---|---|---|---|---|---|---|
| A1 | **CZ CELLxGENE Census LTS 2025-11-08** — chanzuckerberg.github.io/cellxgene-census | Human, mouse, rhesus macaque, marmoset, chimpanzee | scRNA/snRNA (some spatial included, unverified) | 1,845 datasets; human 162,025,130 total / **99,633,637 unique**; mouse 46,299,127 / 21,029,771 unique; macaque 7,010,229; marmoset 2,275,451; chimpanzee 158,099 (Census release page) | Aggregate counts of tissues/cell types/donors not read in this survey (unverified); schema change: the `disease` field can be multi-valued, separated by `' \|\| '` | None (purely observational); some datasets carry disease/developmental-stage metadata | TileDB-SOMA, public S3 bucket; Python `cellxgene_census` ≥1.17 | Datasets mostly CC BY 4.0 (not verified per dataset) | ✅ `gget cellxgene` |
| A2 | **Arc scBaseCount** (Arc Virtual Cell Atlas) — github.com/ArcInstitute/arc-virtual-cell-atlas | 27 species | scRNA (all 10x libraries in SRA, uniformly reprocessed with STARsolo) | **>502M cells**, 27 organisms, 75 tissues (README); older figure on Arc website: 500M / 21 organisms / 72 tissues; published in Cell 2026 | Sample-level metadata includes organism / tissue / disease / **perturbation** / **cell line** plus ontology IDs; SRAgent processed 208,939 SRX, of which 105,343 are 10x (bioRxiv v3) | Includes Velocyto (spliced/unspliced) count features → enables RNA-velocity-style temporal direction; perturbation metadata field (coarse-grained) | `gs://arc-institute-virtual-cell-atlas` (Google Cloud Marketplace; requester-pays for projects that have not subscribed); the old bucket `gs://arc-scbasecount` was shut down on 2026-03-31 | CC0 1.0 (Arc website) | ❌ |
| A3 | **Human Cell Atlas Data Portal** — data.humancellatlas.org | Human (mainly) | Multi-omics, open and managed access | **70.9M cells, 11.3k donors, 532 projects, 1.0k labs** (portal home page statistics) | 2025-12 additions include retina multi-omics, the Asian Immune Diversity Atlas, etc. | A few projects include time/development; no systematic direction | Azul API; AWS Open Data | Varies by project; includes managed access | ❌ (some projects are mirrored into CELLxGENE and reachable indirectly) |
| A4 | **Tabula Sapiens 2.0** — registry.opendata.aws/tabula-sapiens | Human | scRNA (10x + Smart-seq2) | >1.1M cells | **28 tissues, 24 donors (11 M / 13 F, ages 22–74), 182 cell types** (bioRxiv / Cell 2026) | None | AWS Open Data; CELLxGENE | Unverified (CZ-hosted data is generally CC BY 4.0) | ✅ via cellxgene |
| A5 | **Tabula Muris / Muris Senis** | Mouse | scRNA | Unverified | Multiple organs × multiple ages (exact numbers unverified) | Age series (not causal) | CELLxGENE / AWS | Unverified | ✅ via cellxgene (if in Census) |
| A6 | **DISCO** — disco.bii.a-star.edu.sg | Human | Integrated scRNA/multi-omics | >18M cells, 4,593 samples, 351 projects (NAR 2022 version) | **107 tissues/cell lines/organoids, 158 diseases, 20 platforms**; 27 sub-atlases (NAR 2022); an updated NAR 2025 version exists (numbers not read, unverified) | None | `DISCOtoolkit` Python package for metadata-filtered download | Unverified | ❌ |
| A7 | **ARCHS4** (uniformly realigned bulk RNA-seq) — maayanlab.cloud/archs4 | Human, mouse | Bulk | Sample count unverified | All bulk data in GEO; extremely broad contexts | Contains many GEO treatment/time-series experiments (requires metadata mining) | H5 download | Unverified | ⚠️ `gget archs4` provides only correlation/tissue summaries |
| A8 | **CZI Billion Cells Project** (launched 2025-02, with 10x + Ultima) | Human and others | scRNA + perturbation | Target 1B cells, including **250M CRISPR-perturbed cells** (press release) | Plans to cover three categories: natural genetic variation, induced genetic perturbation, chemical perturbation | Yes (planned) | Release status unverified; expected via the CZI Virtual Cells Platform | Stated to be open | ❌ |

---

## B. Perturbation data collections (the core source of within-cell direction)

| # | Name / URL | Organism | Modality | Scale | Diversity | Direction information | Access / size | License | gget |
|---|---|---|---|---|---|---|---|---|---|
| B1 | **scPerturb** — scperturb.org | Human, mouse | scRNA / scATAC / CITE-seq perturbation | 44 datasets (32 CRISPR-RNA, 9 drug-RNA, 3 scATAC, 3 CITE-seq) | >8,000 CRISPR perturbations shared across multiple datasets (paper) | ✅ Intervention (gene/drug) | Zenodo / Figshare+ h5ad (uniform QC + annotation) | Unverified (Zenodo is mostly CC BY 4.0) | ❌ |
| B2 | **PerturBase** — github.com/bm2-lab/PerturBase | Human, mouse | scRNA/multimodal perturbation | 122 datasets / 46 studies, ≈5M cells | **24,254 genetic + 230 chemical perturbations**; 115 single-modal + 7 multimodal (NAR 2025) | ✅ | Web + Zenodo (record 17972932) | Unverified | ❌ |
| B3 | **PerturbDB** (NAR 2025 D1120) | Unverified | Perturbation | Unverified | Unverified | ✅ | Unverified | Unverified | ❌ |
| B4 | **Replogle 2022 GWPS** — gwps.wi.mit.edu; Figshare+ | Human | CRISPRi Perturb-seq | **>2.5M cells** (abstract) | K562 genome-wide 9,866 genes; K562 essential 2,057; RPE1 essential 2,393 (2 cell lines) | ✅ Strong (single-gene knockdown → whole transcriptome) | Figshare+ h5ad; also on CZI VCP and in pertpy | Unverified | ❌ |
| B5 | **Nadig 2025 (TRADE) HepG2/Jurkat essential** | Human | CRISPRi Perturb-seq | Unverified | 2,393 essential-gene perturbations each × 2 cell lines (Nat Genet 2025) | ✅ | Unverified (GEO/Figshare) | Unverified | ❌ |
| B6 | **X-Atlas/Orion** (Xaira) — hf.co/datasets/Xaira-Therapeutics/X-Atlas-Orion | Human | FiCS Perturb-seq (dual-guide CRISPRi) | **~8M cells**; HF API check: 340 entries / **126.3 GB** parquet | 2 cell lines (HCT116, HEK293T) × 18,903 protein-coding genes; median 16,000 UMI/cell; median knockdown efficiency 75.4% / 51.5% (includes dosage information) | ✅ Strong + dosage | HF `datasets` streaming; processed h5ad on Figshare+ (doi 10.25452/figshare.plus.29190726) | **CC BY-NC-SA 4.0** | ❌ |
| B7 | **X-Atlas/Pisces** (Xaira) — hf.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces | Human | CRISPRi Perturb-seq | Claimed **25.6M cells** | Dataset card lists 16 contexts (HCT116, HEK293T, HepG2, iPSC, resting / CD3-CD28-activated Jurkat, multi-lineage differentiated iPSC); press release says 7 contexts (inconsistent figures) | ✅ Strong | **As of 2026-10-04 the HF repo contains only `.gitattributes` / `LICENSE.md` / `README.md` (24 kB); card marked "Coming Soon"; lastModified 2026-03-25.** Planned release: all HCT116/HEK293T; for HepG2/iPSC only a 200-gene validation set + NTC; all resting/activated Jurkat | **CC BY-NC-SA 4.0** | ❌ |
| B8 | **Tahoe-100M** — hf.co/datasets/tahoebio/Tahoe-100M | Human | Drug-perturbation scRNA (Mosaic multi-cell-line pooling) | >100M cells; HF card: 429 GB / 4.3B rows (including DE tables); checked: 3,388 parquet shards under data/ | **50 cancer cell lines × ~1,100 small molecules**, multiple doses; metadata: 1.34k samples, 379 drugs (with targets/MOA), 50 lines (with driver mutations), pseudobulk DE 4.09B rows | ✅ (chemical intervention; targets indirect) | HF parquet streaming; also in Arc VCA | **CC0 1.0** | ❌ |
| B9 | **Parse 10M PBMC cytokines / Human Cytokine Dictionary** — parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment | Human | Ligand (cytokine) stimulation scRNA | **9,697,974 cells** | **12 donors × (90 cytokines + PBS) = 1,092 conditions**; 18 immune cell types; 24 h; median 7,400 cells/condition | ✅ **Between-cell direction** (ligand → response of each cell type) | Registration form required; h5ad 41 GB, pseudobulk 6 GB, DEG 3.2 GB; HF mirror `slaf-project/Parse-10M` (SLAF format) | **CC BY-NC 4.0** | ❌ |
| B10 | **Primary human CD4+ T cell genome-scale Perturb-seq** (Zhu, Dann, …, Pritchard, Marson; Cell 2026) — virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq | Human (primary) | CRISPRi Perturb-seq (probe-based) | **~22M cells** | **4 donors × 3 conditions (Rest / Stim 8 h / Stim 48 h)**; all expressed genes; DE tables cover 33,983 perturbations; DE stratified by guide/donor (h5mu) | ✅ Strong + **stimulation time dimension** + donor differences | S3 (AWS CLI / VCP CLI); cell-level `D*_*.assigned_guide.h5ad`, `GWCD4i.pseudobulk_merged.h5ad` | **MIT** | ❌ |
| B11 | **Jiang, Dalgarno, Papalexi, Satija 2025 (Nat Cell Biol)** — GEO GSE281048; Zenodo 14518762 | Human | Perturb-seq × signaling stimulation | Unverified | **>1,500 perturbations × 6 cell lines (A549, MCF7, HT29, HAP1, BxPC3, K562) × 5 stimuli (IFNB, IFNG, TGFB, TNFA, INS)** | ✅ **Signaling-pathway direction** (receptor/regulator knockout × ligand stimulation) | GEO + Zenodo (Seurat objects, DE results) | Unverified | ❌ |
| B12 | **LINCS L1000 / CMap 2020** — clue.io/releases/data-dashboard; GEO GSE92742/GSE70138 (accessions unverified) | Human | Bulk L1000 (978 landmark genes + inferred) | Level 5: **921,123 signatures** | **248 cell lines** (L2S2 figure); 12,487 compounds, 4,345 shRNA, 5,157 CRISPR, 4,040 overexpression perturbations; 10,174 high-quality genes in total | ✅ (genetic + chemical intervention; broadest across cell lines) | GCTX download; processed versions on Zenodo | Unverified | ⚠️ Possibly only via Enrichr's LINCS gene-set libraries for enrichment (library names unverified) |
| B13 | **DepMap / CCLE** — depmap.org/portal/data_page | Human (cancer cell lines) | Bulk RNA/WGS + CRISPR KO dependency | Total model count not read in this survey (unverified; ~1,000+ is common knowledge, unverified) | 25Q2 +5 screens, 25Q3 +3 screens, **26Q1 +25 genome-wide screens** (release notes); no longer posted to Figshare after 24Q4 | ⚠️ Dependency (gene → fitness), not gene → gene | Manual / API download from the portal | DepMap terms (unverified) | ⚠️ `gget opentargets resource='depmap'` gives only gene-level summaries |
| B14 | **Open Problems OP3 (NeurIPS 2023 perturbation prediction)** — openproblems.bio/benchmarks/perturbation_prediction; GEO GSE279945 | Human (primary PBMC) | Drug-perturbation scRNA (+ partial multiome) | Unverified | **146 compounds × 3 donors**, 24 h, multiple immune cell types | ✅ (chemical; across cell types) | GEO / S3 (unverified) | Unverified | ❌ |

---

## C. Time series / development / differentiation / lineage / RNA kinetics (temporal direction)

| # | Name / URL | Organism | Modality | Scale | Diversity | Direction information | Access | License | gget |
|---|---|---|---|---|---|---|---|---|---|
| C1 | **Qiu et al. 2024 Nature, mouse embryo time-lapse** — CELLxGENE collection 45d5d2c3… | Mouse | sci-RNA-seq3 snRNA | 12.4M nuclei / 83 embryos sequenced; after filtering **11,441,407 nuclei / 74 embryos** | **E8 → P0, at 2 h (somite stages) to 6 h intervals; 20-min resolution right after birth** | ✅ Time series (development) | CELLxGENE (author website also exists, unverified) | Unverified | ✅ via cellxgene (if in Census, unverified) |
| C2 | **Cao 2020 human fetal gene expression atlas** | Human | sci-RNA-seq3 | **4,062,980 cells, 121 samples** | **15 organs, 10–18 weeks post-fertilization** | ⚠️ Weak time (gestational age) | Author website / GEO (unverified) | Unverified | Unverified |
| C3 | **Cuomo 2020 iPSC → endoderm** | Human | Smart-seq2 scRNA | 36,044 cells | **125 iPSC lines / 125 donors × 4 time points**; >30% of eQTLs specific to a single stage | ✅ Time + natural genetic variation (dynamic eQTL) | Unverified (ENA/Zenodo) | Unverified | ❌ |
| C4 | **Jerber 2021 iPSC → dopaminergic neurons** | Human | scRNA | **>1M cells** | **215 iPSC lines × 3 differentiation time points** + rotenone oxidative stress | ✅ Time + genetics + stimulus | Unverified | Unverified | ❌ |
| C5 | **scLTdb (single-cell lineage tracing database)** — scltdb.com | Human, mouse, etc. (3 species) | Lineage-tracing scRNA/scATAC | **109 datasets, 2.8M cells** | **36 scLT technologies** (LARRY, CellTag, SISBAR, …), 13 categories of tissue origin (HSPC 37, cell lines 22, brain 22 datasets) | ✅ Clonal lineage (state → fate) | h5ad / Seurat download | Unverified | ❌ |
| C6 | **LARRY (Weinreb 2020 Science)** | Mouse HSPC | Lineage-barcoded scRNA | Unverified | In vitro / in vivo, multiple time points; **2,632 clones across multiple time points** | ✅ Lineage + time | scLTdb / GEO | Unverified | ❌ |
| C7 | **CellTag-multi (Jindal 2023 NBT)** — GEO GSE216506 | Mouse | Lineage + scRNA + scATAC | Unverified | Hematopoiesis + fibroblast → induced endoderm progenitor reprogramming | ✅ Lineage + chromatin → expression | GEO | Unverified | ❌ |
| C8 | **Metabolic-labeling scRNA (RNA kinetics)**: sci-fate (A549 + dexamethasone 0/2/4/6/8/10 h), scNT-seq (mouse cortical neurons, KCl 0/15/30/60/120 min), scEU-seq (GSE128365, organoids + cell cycle) | Human/mouse | 4sU nascent-RNA labeling | Small (each roughly thousands to ~100k cells, unverified) | Stimulation time series; synthesis/degradation rates | ✅ **New vs. old RNA → arrow of time** (stronger than splicing velocity) | GEO; built-in downloaders in the `dynamo` package | Unverified | ❌ |
| C9 | **PerturbSci-Kinetics** (Nat Biotech 2023) | Human | CRISPR perturbation × metabolic labeling | Unverified | Unverified | ✅✅ Intervention + kinetics | Unverified | Unverified | ❌ |
| C10 | **CD4 T stimulation time course (Rest / 8 h / 48 h from B10)** | — | — | — | See B10 | ✅ Perturbation × time | — | MIT | ❌ |

---

## D. Spatial and cell–cell communication (between-cell direction)

| # | Name / URL | Organism | Modality | Scale | Diversity | Direction information | Access | License | gget |
|---|---|---|---|---|---|---|---|---|---|
| D1 | **SpatialCorpus-110M (Nicheformer, Nat Methods 2025)** | Human, mouse | Dissociated scRNA + imaging-based spatial (MERFISH, Xenium, CosMx, ISS) | **>110M cells = 57.06M dissociated + 53.83M spatial** | Dissociated part: 73 tissues; spatial part: **15 solid organs**; includes 18 cell lines (from search-result abstract) | ⚠️ Spatial neighborhood (correlational; needs LR priors to obtain direction) | The paper provides data-collection code at `github.com/theislab/nicheformer-data` (returned 404 when accessed); **no packaged corpus download found (unverified)** | Each source dataset keeps its own license | ❌ |
| D2 | **HEST-1k** (Mahmood Lab, NeurIPS 2024) — github.com/mahmoodlab/HEST; HF | Human, mouse | ST + H&E whole-slide images | **1,229 samples, 2.1M spots, >76M cells (nucleus segmentation)** | 153 cohorts, **26 organs, 25 cancer types (367 cancer samples)**; Visium 49.0%, STv1 44.9%, Xenium 5.3%, Visium HD 0.8% | ⚠️ Spatial neighborhood | HF download (may require accepting terms, unverified) | Unverified | ❌ |
| D3 | **Allen Brain Cell (ABC) Atlas — whole mouse brain MERFISH** — alleninstitute.github.io/abc_atlas_access | Mouse | MERFISH (500 genes) + sc/snRNA | **~4M spatial cells**, registered to CCFv3 | **34 classes / 338 subclasses / 1,201 supertypes / 5,322 clusters** | ⚠️ Spatial (+ a separate MERFISH dataset from the Zhuang lab) | AWS S3 + `abc_atlas_access` tutorials | Allen terms (unverified) | ❌ |
| D4 | **HTAN (including the HTAPP pilot)** — humantumoratlas.org | Human | scRNA, scATAC, multiplexed imaging (MIBI/CODEX/MxIF), spatial | **2,464 participants, 10,626 biospecimens** (Release 6.0); ~30,000 files (Release 4.0) | 10 centers, multiple cancer types + precancerous lesions | ⚠️ Tumor-microenvironment spatial neighborhood | Open: HTAN Portal, ISB-CGC BigQuery; controlled: CRDC/DRS | Open + controlled | ❌ |
| D5 | **Immune Dictionary (Cui 2024 Nature)** — Broad SCP2554 | Mouse (in vivo lymph node) | In vivo cytokine-injection scRNA | Unverified | **86 cytokines × >17 immune cell types (>1,400 combinations)**; >66 cytokine-driven polarization states | ✅ **In vivo ligand → cell-type response** | Single Cell Portal | Unverified | ❌ |
| D6 | **Spatial perturbation** (Spatial Perturb-seq / Stereo-seq, Nat Commun 2026; Perturb-DBiT, GSE319277; PerturbSpace, bioRxiv 2026) | Mostly mouse | In vivo CRISPR + spatial | Small | Mouse brain neurodegeneration risk genes; 40 transcriptional regulators in splenic hematopoiesis; liver immune niches | ✅✅ **Intervention + neighborhood → non-cell-autonomous effects** | GEO | Unverified | ❌ |
| D7 | **CellPhoneDB v5** — github.com/ventolab/cellphonedb-data | Human | LR prior (including multi-subunit complexes) | **~2,912 interactions**, ~1,000 of them with non-peptide ligands (hormones, neurotransmitters, etc.) | — | ✅ Ligand → receptor direction | GitHub | Unverified | ❌ |
| D8 | **CellChatDB v2** | Human, mouse | LR prior | **~3,300 interactions** (~40% secreted, 17% ECM–receptor, 13% contact, ~30% non-protein) | — | ✅ | CellChat R package | Unverified | ❌ |
| D9 | **OmniPath intercell + LIANA(+)** — omnipathdb.org | Human (orthology-mapped to mouse) | Integrated LR / intercell role annotation | OmniPath curated LR ~6,307 (per the OmnipathR function) | Role classes: ligand / receptor / transporter, etc. | ✅ | REST API; `liana`, `omnipath` packages | Annotated per source resource; most allow commercial use (NAR 2026) | ❌ |
| D10 | **NicheNet v2 prior** (ligand → receptor → signaling → TF → target gene) | Human, mouse | Integrated network | Roughly twice as many ligands as v1 (exact number unverified) | — | ✅ Ligand → **downstream target genes** (a direction chain from between-cell to within-cell) | Zenodo / nichenetr | Unverified | ❌ |
| D11 | **Lignature** (Genome Res 2026) / **CytoSig** (Nat Methods 2021) | Human | Ligand-induced transcriptomic signatures (collected from public experiments) | Lignature: **362 ligands**; CytoSig: 20,591 cytokine/chemokine/growth-factor response experiments, reliably predicts the activity of 43 cytokines | — | ✅ Ligand → response signature | Website / GitHub | Unverified | ❌ |

---

## E. Directed gene–gene knowledge sources (weak supervision / evaluation sets)

| # | Name / URL | Scale | Direction information | Access | License | gget |
|---|---|---|---|---|---|---|
| E1 | **OmniPath (NAR 2026)** — omnipathdb.org | Integrates **168 resources**; interactions domain **1,419,006 unique interactions / 115 resources** (signaling, TF–target, miRNA, drug–target); enzyme–PTM 115,215; complexes 52,086 | ✅ Directed + signed (stimulation/inhibition), with a majority-vote "consensus" column | REST API (JSON/TSV) | Annotated per resource; "most allow commercial use" | ❌ |
| E2 | **SIGNOR 4.0 (NAR 2026)** | >8,800 new entries (+27%); total not read in this survey; extrapolating from ~33,000 entries in 3.0 gives ~42,000 (**extrapolated, unverified**) | ✅ Manually curated causal relations (directed, signed) | Website download / API (details unverified) | Unverified | ❌ (possibly partly visible via opentargets interactions, unverified) |
| E3 | **CollecTRI** (NAR 2023) | **1,183 TFs, 45,856 signed TF–gene interactions** | ✅ TF → target gene (signed) | Via OmniPath / decoupler | Unverified | ❌ |
| E4 | **DoRothEA** | Human: **1,395 TFs → 20,244 genes, 486,676 interactions**; A–E confidence levels | ✅ TF → target gene | Bioconductor / OmniPath | Unverified | ❌ (gget enrichr's ChEA_2016 is only ChIP-derived gene sets) |
| E5 | **Reactome V96** — reactome.org | **16,338 human reactions, 2,870 pathways, 32,399 proteins (11,452 genes), 16,145 complexes** | ✅ Reaction-level direction (input → output, catalysis, regulation) | REST / Neo4j graph database download | Unverified (officially usually CC BY 4.0) | ❌ (gget PR #236 was closed; enrichment possible via Enrichr's Reactome library) |
| E6 | **ENCODE-rE2G / enhancer–gene encyclopedia** | **13.5M enhancer–gene regulatory pairs across 352 cell types/tissues**; benchmark of **10,411 CRISPR-tested element–gene pairs** (plus >30,000 fine-mapped eQTLs and 569 GWAS variants) | ✅ Cis direction (enhancer → gene), CRISPRi gold standard | ENCODE portal | ENCODE open | ❌ (gget PR #231 encode was closed) |
| E7 | **GTEx v10** — gtexportal.org | **943 donors, 19,466 samples (eQTL), 50 tissues**; 16,760 new smallRNA-seq samples | ⚠️ Natural genetic variation → cis-eQTL (MR-style direction anchor) | Portal API v2 / AnVIL | Open tier + dbGaP controlled | ⚠️ `gget opentargets` has only baseline expression; no eQTLs |
| E8 | **OneK1K (single-cell eQTL)** — GEO GSE196830; HCA | **1.27M PBMCs, 982 donors, 14 immune cell types**; 26,597 cis-eQTLs, 990 trans-eQTLs | ✅ trans-eQTL = "genetically anchored gene → gene" direction (cell-type-specific) | GEO / HCA | Unverified | ❌ |
| E9 | **inspre causal network (Brown et al., Nat Commun 2025)** | **788 genes** from the Replogle K562 GWPS; **131,155 significant effects** at 5% FDR | ✅ Directed network inferred from interventional data (small-world, scale-free) | Paper supplement / GitHub (unverified) | Unverified | ❌ |

---

## Synthesis

### (1) Sources with the highest context diversity per byte

Rough ranking by "number of independent contexts (cell type × tissue × donor × condition × cell line × time) / download size" (sizes are mostly order-of-magnitude estimates; unverified points are marked):

1. **LINCS L1000 Level 5**: 248 cell lines × ~25,000 perturbations, ~920k signatures, only 978 measured genes, very small size → the highest number of "cell line × perturbation" combinations per byte. Drawbacks: bulk, narrow gene coverage.
2. **Pseudobulk / DE products of perturbation datasets**: CD4 T `GWCD4i.pseudobulk_merged.h5ad` (4 donors × 3 conditions × genome-wide perturbations), Tahoe pseudobulk DE tables, Parse pseudobulk (6 GB, covering 1,092 conditions). These compress "number of cells" into "number of conditions", which fits "diversity over quantity" exactly.
3. **CELLxGENE Census (unique cells + stratified downsampling by dataset × cell type × tissue)**: 1,845 datasets across 5 species. Because SOMA supports slicing by metadata, one can pull only up to a cap of cells per context.
4. **scBaseCount**: the broadest contexts (27 species / 75 tissues, with perturbation and cell-line metadata), but also the largest volume (>500M cells); requires context deduplication on sample metadata before downsampling; request fees apply on GCS.
5. **GTEx v10 / ARCHS4 bulk**: 50 tissues × 943 donors / all of GEO; small size, high donor and tissue diversity.
6. **Tabula Sapiens 2.0**: the same donors across 28 tissues, 182 cell types; clean and controlled, but only ~1.1M cells.

Inefficient sources (low diversity per byte): **very large perturbation sets from a single cell line under a single condition** (e.g., using only the HCT116 portion of Orion; 126 GB covers only 2 cell lines). Such data should be sampled per condition rather than fed in full.

### (2) Ranking by direction information (strongest to weakest)

1. **Multi-context genome-scale genetic perturbation (within-cell gene → gene)**: CD4 T (primary, 4 donors × 3 stimulation time points) > X-Atlas/Orion (with dosage) > Replogle K562/RPE1 > Nadig HepG2/Jurkat > X-Atlas/Pisces (not yet released).
2. **Perturbation × signaling stimulation (pathway direction)**: Jiang 2025 (6 cell lines × 5 ligands × >1,500 perturbations). This is the only public collection that systematically crosses "receptor/regulator knockout" with "ligand stimulation", carrying both within-cell and between-cell direction at once.
3. **Ligand as intervention (between-cell direction)**: Parse 10M (90 cytokines × 12 donors × 18 cell types), Immune Dictionary (in vivo, 86 cytokines), Lignature / CytoSig (signature level).
4. **Spatial perturbation** (Spatial Perturb-seq, Perturb-DBiT, PerturbSpace): intervention plus neighborhood, able to separate non-cell-autonomous effects; the most direct direction signal, but small scale.
5. **Intervention + kinetics**: PerturbSci-Kinetics; **metabolic-labeling time series** (sci-fate, scNT-seq, scEU-seq): nascent RNA provides an arrow of time.
6. **Chemical perturbation**: Tahoe-100M, OP3, sci-Plex (via scPerturb), L1000 compounds. These are interventions, but targets are indirect; they must be combined with drug–target priors.
7. **Natural genetic variation**: OneK1K trans-eQTLs, Cuomo/Jerber dynamic eQTLs, GTEx cis-eQTLs. Usable for MR-style direction anchoring.
8. **Lineage / developmental time**: scLTdb, LARRY, CellTag-multi, Qiu 2024. These give state → fate and temporal order; weakly causal.
9. **Spatial neighborhood + LR priors** (SpatialCorpus, HEST, ABC, HTAN): correlational data; direction comes from the LR priors.
10. **Knowledge bases** (OmniPath, SIGNOR, CollecTRI, DoRothEA, Reactome, ENCODE-rE2G): directed labels, suited for weak supervision or evaluation, not as the main training signal; among them, the 10,411 ENCODE CRISPRi pairs are the gold standard for cis direction.

### (3) gget coverage gaps and the most valuable new gget modules

**Current state**: of the ~45 sources above, the only one gget can access directly is CELLxGENE Census (including subsets such as Tabula Sapiens and Qiu 2024). ARCHS4, DepMap, PPIs, and LINCS / ChEA / Reactome are reachable only as summary statistics or enrichment results via `archs4`, `opentargets`, and `enrichr`. **None of the perturbation, spatial, lineage, eQTL, or directed signaling-network resources are in gget.**

Proposed new modules (ranked by "value to this project × fit with gget's data-access role × stable API that allows live tests"):

| Priority | Proposed module | Rationale | Feasibility |
|---|---|---|---|
| 1 | **`gget omnipath`** (interactions / intercell / enz_sub / complexes; optional CollecTRI, DoRothEA, SIGNOR subsets) | A single REST endpoint yields directed, signed signaling networks, TF regulation, and LR pairs, filling exactly the gap in directional priors; results are lightweight tables, consistent with gget's style | High: REST, no authentication, easy to write TDD live tests (e.g., upstream regulators of `TP53`) |
| 2 | **`gget perturb`** (perturbation data catalog and download: scPerturb's Zenodo manifest, CZI Virtual Cells Platform datasets, Tahoe / Orion / Parse mirrors on HF; query by cell line, perturbation type, target gene; return URLs and metadata, optionally stream subsets) | Perturbation data is currently scattered across five places (Zenodo, Figshare, S3, HF, GEO), which is the biggest access pain point | Medium: requires maintaining a metadata index; HF and VCP have APIs |
| 3 | **`gget scbasecount`** (query sample-level metadata by organism / tissue / perturbation / cell_line, return GCS paths) | The corpus with the broadest contexts, and it carries perturbation and cell-line fields | Medium: GCS requester-pays billing must be flagged in the docs |
| 4 | **`gget gtex`** (GTEx API v2: eQTL, sQTL, tissue expression) | Genetically anchored direction information; stable API | High |
| 5 | **`gget lincs`** (SigCom LINCS / L2S2 API: search signatures by gene or perturbation) | The perturbation resource with the highest diversity per byte | Medium-high: Ma'ayan Lab API, same origin as enrichr/archs4 (see gget issue #164) |
| 6 | Revive `gget reactome` (PR #236) and `gget encode` (PR #231, focusing on CRISPRi / rE2G files) | Prototype code already exists | First need to understand why they were closed (unverified) |
| 7 | `gget depmap` (issue #121 closed) | Partly superseded by opentargets; low incremental value | Low priority |

### (4) Phased pretraining data plan

**Phase 0: vocabulary and priors (1–2 weeks)**
- Unify gene IDs (Ensembl, with human–mouse orthology mapping). Build a directed prior graph: OmniPath / CollecTRI / SIGNOR / CellPhoneDB / NicheNet.
- These are for evaluation and optional auxiliary losses, **not the bulk of training labels**. Also define held-out sets at this stage: held-out cell lines, donors, perturbed genes, and ligands.

**Phase 1: observational diversity pretraining (sampling by context)**
- Core: Census LTS 2025-11-08 (unique cells) plus the non-overlapping part of scBaseCount.
- Define the sampling unit as `(dataset, tissue, cell_type, donor/cell_line, disease, stage)`, with **a cell cap per unit** (e.g., a few hundred). This increases the number of contexts rather than the number of cells.
- Add bulk contexts: GTEx (50 tissues × 943 donors) and ARCHS4, as pseudo-cell tokens.
- Add spatial contexts: the sources in SpatialCorpus, HEST-1k, ABC MERFISH, the HTAN open tier. Introduce "neighborhood tokens" during training to pave the way for between-cell modeling.

**Phase 2: mid-training with directional signals**
- Time and kinetics: Qiu 2024 (E8–P0), Cao 2020 fetal atlas, iPSC differentiation (Cuomo, Jerber), metabolic labeling (sci-fate / scNT / scEU). scBaseCount's Velocyto layers can supply unspliced/spliced information. This part trains a temporal prediction objective, i.e., predicting t+Δ from t.
- Lineage: the 109 datasets in scLTdb, training a state → fate objective.
- Genetic anchors: eQTLs from OneK1K and GTEx, as an auxiliary task on "variant → gene → downstream" triples.
- Between-cell: spatial neighborhoods plus LR priors, training masked prediction of "sender expression → receiver state".

**Phase 3: perturbation fine-tuning (interventional supervision, evaluated on held-out contexts)**
- Genetic perturbation: CD4 T (4 donors × 3 conditions), Orion (2 cell lines, with dosage), Replogle (K562/RPE1), Nadig (HepG2/Jurkat). Add X-Atlas/Pisces once it is public.
- Signaling × genetics: Jiang 2025 (6 cell lines × 5 ligands).
- Between-cell intervention: Parse 10M, Immune Dictionary; finally add spatial perturbation data (Perturb-DBiT, etc.).
- Chemical perturbation: Tahoe-100M (50 cell lines), OP3, L1000 (248 cell lines). Drugs must be tokenized using target priors.
- Evaluation: **hold out by unseen cell line, donor, and stimulation condition.** Use the category E knowledge bases and ENCODE CRISPRi to check predicted directions.

**License note**: Orion and Pisces are CC BY-NC-SA 4.0 and Parse is CC BY-NC 4.0; all are non-commercial licenses. Tahoe and scBaseCount are CC0. CD4 T is MIT. If commercialization is intended, non-commercially licensed data should be isolated in a separate training branch.

---

## Sources

- Census data releases: https://chanzuckerberg.github.io/cellxgene-census/cellxgene_census_docsite_data_release_info.html
- scBaseCount README: https://github.com/ArcInstitute/arc-virtual-cell-atlas/blob/main/scBaseCount/README.md ; bioRxiv: https://www.biorxiv.org/content/10.1101/2025.02.27.640494v3 ; Arc VCA: https://arcinstitute.org/tools/virtualcellatlas
- HCA Data Portal: https://data.humancellatlas.org/
- Tabula Sapiens 2.0: https://www.cell.com/cell/fulltext/S0092-8674(26)00937-2 ; https://www.biorxiv.org/content/10.1101/2024.12.03.626516v2.full
- DISCO: https://academic.oup.com/nar/article/50/D1/D596/6430491 ; https://academic.oup.com/nar/article/53/D1/D932/7899529
- CZI Billion Cells: https://www.biospace.com/press-releases/chan-zuckerberg-initiative-launches-billion-cells-project-with-10x-genomics-and-ultima-genomics-to-advance-ai-in-biology
- scPerturb: https://www.biorxiv.org/content/10.1101/2022.08.20.504663v3 ; https://plus.figshare.com/articles/dataset/scPerturb_Single-Cell_Perturbation_Data_RNA_and_protein_h5ad_files/24160713
- PerturBase: https://pmc.ncbi.nlm.nih.gov/articles/PMC11701531/ ; PerturbDB: https://academic.oup.com/nar/article/53/D1/D1120/7755478
- Replogle 2022: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9380471/ ; https://gwps.wi.mit.edu/
- Nadig 2025: https://www.nature.com/articles/s41588-025-02169-3
- X-Atlas/Orion: https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion (README + HF API check)
- X-Atlas/Pisces: https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces (HF API check shows only 3 files); https://www.businesswire.com/news/home/20260317710096/en/
- Tahoe-100M: https://huggingface.co/datasets/tahoebio/Tahoe-100M
- Parse 10M PBMC: https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/ ; Human Cytokine Dictionary: https://www.biorxiv.org/content/10.64898/2025.12.12.693897v1
- CD4 T Perturb-seq: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq ; https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
- Jiang 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC12083445/ ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281048
- LINCS L1000: https://pmc.ncbi.nlm.nih.gov/articles/PMC12230732/ ; https://clue.io/releases/data-dashboard
- DepMap release notes: https://forum.depmap.org/t/depmap-quarterly-release-notes/3560
- OP3: https://openproblems.bio/benchmarks/perturbation_prediction/ ; https://proceedings.neurips.cc/paper_files/paper/2024/file/24c4d51f3ef48dd2dbab78243ecb26a1-Paper-Datasets_and_Benchmarks_Track.pdf
- Qiu 2024: https://www.nature.com/articles/s41586-024-07069-w
- Cao 2020 fetal atlas: https://www.science.org/doi/10.1126/science.aba7721
- Cuomo 2020: https://www.nature.com/articles/s41467-020-14457-z ; Jerber 2021: https://www.nature.com/articles/s41588-021-00801-6
- scLTdb: https://academic.oup.com/nar/article/53/D1/D1173/7848853 ; LARRY: https://www.science.org/doi/10.1126/science.aaw3381 ; CellTag-multi: https://www.nature.com/articles/s41587-023-01931-4
- sci-fate: https://pmc.ncbi.nlm.nih.gov/articles/PMC7416490/ ; scNT-seq: https://pmc.ncbi.nlm.nih.gov/articles/PMC8103797/ ; dynamo: https://www.cell.com/cell/fulltext/S0092-8674(21)01577-4 ; PerturbSci-Kinetics: https://www.nature.com/articles/s41587-023-01948-9
- Nicheformer / SpatialCorpus-110M: https://pmc.ncbi.nlm.nih.gov/articles/PMC12695652/
- HEST-1k: https://arxiv.org/abs/2406.16192 ; https://github.com/mahmoodlab/hest
- ABC Atlas: https://alleninstitute.github.io/abc_atlas_access/ ; https://www.nature.com/articles/s41586-023-06812-z
- HTAN: https://pmc.ncbi.nlm.nih.gov/articles/PMC12668157/
- Immune Dictionary: https://www.nature.com/articles/s41586-023-06816-9
- Spatial Perturb-seq: https://www.nature.com/articles/s41467-026-69677-6 ; Perturb-DBiT: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12083649/ ; PerturbSpace: https://www.biorxiv.org/content/10.64898/2026.05.25.727765v1.full
- CellPhoneDB v5: https://arxiv.org/pdf/2311.04567 ; CellChat v2: https://www.biorxiv.org/content/10.1101/2023.11.05.565674v1.full
- OmniPath NAR 2026: https://academic.oup.com/nar/article/54/D1/D652/8326458 ; https://r.omnipathdb.org/reference/curated_ligand_receptor_interactions.html
- NicheNet: https://www.nature.com/articles/s41596-024-01121-9 ; Lignature: https://pmc.ncbi.nlm.nih.gov/articles/PMC13445700/ ; CytoSig: https://www.nature.com/articles/s41592-021-01274-5
- SIGNOR 4.0: https://academic.oup.com/nar/article/54/D1/D682/8324960
- CollecTRI: https://academic.oup.com/nar/article/51/20/10934/7318114 ; DoRothEA: https://saezlab.github.io/dorothea/articles/dorothea.html
- Reactome V96: https://reactome.org/about/news/291-v96-released
- ENCODE-rE2G: https://pmc.ncbi.nlm.nih.gov/articles/PMC10680627/
- GTEx v10: https://anvilproject.org/news/2024/11/20/gtexv10
- OneK1K: https://www.science.org/doi/10.1126/science.abf3041
- inspre: https://www.nature.com/articles/s41467-025-64353-7
- gget source code and issues/PRs: https://github.com/scverse/gget (`gget/` directory, issues #218/#121/#164, PRs #231/#236 checked)
