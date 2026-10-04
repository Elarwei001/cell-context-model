# 2026 Literature Review: Zero-Shot Cross-Cell-Line Genetic Perturbation Prediction × In-Context Learning / Cross-Cell Attention × Flow Matching

**Date:** 2026-10-04

**Scope.** This review surveys work published up to early October 2026 on predicting single-cell transcriptional responses to genetic perturbations in cellular contexts that were never seen with any perturbation during training. The motivating problem is the core evaluation question of this project, which studies how genes and cells influence each other across diverse cellular contexts: predicting perturbation responses in cellular contexts never seen perturbed — e.g. zero-shot CRISPRi prediction in a cell line or cell type for which only unperturbed, non-targeting control (NTC) cells are available. In the regime of most interest, targets are mostly non-essential genes with weak effects, and the goal is to generate realistic single-cell count profiles for each target in the new context. The question driving the review is whether a control-anchored flow-matching approach we explored earlier (conditional flow matching from control to perturbed cells, with a cell-line embedding computed only from that line's control set, a fixed stratified gene axis, and emission on real control counts) is worth extending into an in-context learning (ICL) model: concatenate related perturbed cells (from other cell lines) with the new line's control cells into a single input and use bidirectional cross-cell attention (set / axial / MSA-style row-column attention, implemented with varlen packing or FlexAttention). Literature search date: 2026-10-04.

**Note on source access.** bioRxiv full texts were largely inaccessible during this survey (see Section 0), so several bioRxiv items rely on abstracts and metadata only. These are marked explicitly in the entries below.

---

## 0. Reading notes (level of verification)

- **bioRxiv web pages were blocked by Cloudflare** (both WebFetch and curl returned 403). For bioRxiv papers, only abstracts / authors / dates / licences were read via `api.biorxiv.org`. Where a PMC version exists (Stack v1, Virtual Cells Need Context, pertTF), the full text was read through Europe PMC full-text XML. X-Cell's full text was read from the PDF on Xaira's CDN.
- arXiv papers were read via HTML/PDF (some were summarised by WebFetch's small model; key numbers were cross-checked by grepping the PDF text wherever possible; those not cross-checked are marked "(paraphrased from a summarisation tool)").
- Nature-family papers (TxPert, NBT 2026; TRADE, NG 2025) were behind a login wall; only the arXiv version / Europe PMC abstract was read.
- Anything **unverified** is explicitly marked `(unverified)`. Numbers that were not actually read are never reported.

---

## 1. Most relevant papers at a glance (ordered by relevance to "ICL for zero-shot cross-line genetic perturbation")

| # | Paper | Date | Keywords | Distance from our setting |
|---|---|---|---|---|
| 1 | Stack (Arc) | 2026-01 / v2 2026-06 | Tabular attention (intra-cell + inter-cell), prompt/query cell-set ICL | Closest architecture; but ICL evaluation is mostly cytokines/drugs/donors, and genetic-perturbation performance is very weak |
| 2 | PT-RAG | 2026-03 | Retrieves related perturbations as context on top of STATE, differentiable sparse selection | Closest idea (related perturbations as context); but the target line is few-shot (30% of perturbations visible) |
| 3 | PerturbMap | 2026-07 | Reliability-weighted source→recipient transfer of the same perturbation | Cross-context transfer + "harm control" idea; requires paired anchors in the target context (not zero-shot) |
| 4 | Tabular FMs (TabPFN/TabICL) are competitive | 2026-07 | General-purpose tabular ICL for cross-cell-type perturbation prediction | Shows that "generic ICL" can match specialised models |
| 5 | X-Cell / X-Atlas-Pisces (Xaira) | 2026-03 | 25.6M CRISPRi cells, 16 contexts; zero-shot + test-time adaptation on control cells | Data + TTA can be borrowed directly |
| 6 | Virtual Cells Need Context, Not Just Scale | 2026-02 | Context diversity > cell count; causal transportability | Principle for organising training data |
| 7 | Mechanisms Matter (Microsoft) | 2026-07 | Cross-context generalisation gap; semi-synthetic causal simulator | Negative evidence |
| 8 | CellMSA | 2026-09 | MSA-style context cells → gene-pair representation → attention bias | Direct precedent for MSA / row-column attention; but only tests unseen perturbations |
| 9 | STATE (Arc) | 2025-06 / Cell 2026 | Cell-set transformer, MMD | Our main reference baseline |
| 10 | C3TL | 2026-03 | Perturbation encoder with mean aggregation across contexts | Shows that "simple set pooling" already matches STATE |
| 11 | scPILOT (Adv Sci) | 2026-08 | Query-conditioned latent OT transfer | Cell-line genetic perturbation benchmark ≈ identity |
| 12 | TxPert | 2025-05 / NBT 2026 | Multi-knowledge-graph GNN; includes a cross-line setting where "the test line has no perturbations at all" | KG prior + cross-line zero-shot |
| 13 | PRiMeFlow (Altos) | 2026-04 | FM in raw gene space + U-Net | Reference for generator design |
| 14 | scDFM (ICLR 2026) | 2026-02 | CFM + MMD + differential attention | Reference for generator design |
| 15 | GeneGeoFlow (control-anchored residual FM) | 2026-08 | Control-anchored residual flow + Δ-correlation loss | Isomorphic to the control-anchored flow-matching approach we explored |
| 16 | scBIG (module-inductive) | 2026-02 | Data-driven gene modules + CFM | Relevant to a "co-response module" axis |
| 17 | TRADE (Nadig 2025) | 2025-04 | A typical perturbation affects ~45 genes; moderate cross-cell-type consistency | The reality of weak effects |

(Full entries in Section 2.)

---

## 2. Paper entries

### 2.A Cross-context / in-context / cells-as-tokens / retrieval augmentation

#### A1. Stack: In-Context Learning of Single-Cell Biology
- **Authors/affiliations**: Dong M., Adduri A., Gautam D., … Burke D.P., Roth T., Roohani Y.H. (Arc Institute / Stanford, etc.)
- **Date/links**: v1 2026-01-09, v2 2026-06-08; bioRxiv 10.64898/2026.01.09.698608 (CC-BY preprint); PMC12803207 (v1 full text)
- **Method**: Each input is a cell set from a single sample. Each cell is projected by a single layer into 100 "gene module tokens" (trainable grouped tokenisation, independent of external gene semantics). Each block = intra-cell MHA (over one cell's module tokens) + inter-cell MHA (cell token = concatenation of all module tokens, attention over the cell set) + FFN, inspired by TabPFN/TabICL. Pre-training: rectangular masking (the same set of genes is masked for all cells at once), reconstructing every cell; plus a distributional regulariser (cell embedding = set constant + standard normal, for linear identifiability). Post-training: the same sample is split by cell type into prompt (visible) and target (hidden); the target is replaced with query cells of the same type from **another sample**, trained with self-distillation (EMA teacher) + distribution matching in expression space under a zero-inflated normal, learning "the counterfactual state of the query conditioned on the prompt". At inference, query cells are replaced iteratively in order of confidence (similar to masked diffusion).
- **Evaluation setting**: 4 classes of ICL tasks: (1) perturbation effects across cell types (prompt = perturbed cells of a random type, query = controls of other types); (2) cross-dataset perturbation; (3) donor / observational conditions; (4) cross-dataset generation of cell types. Data: OpenProblems drugs, Dong 2023 cytokines, immune ageing, Tabula Sapiens, etc. (**mostly cytokines/drugs/donors**). Post-training data: CELLxGENE + Parse PBMC 10M (12 donors, 90 cytokines).
- **Key results** (v1 main text): zero-shot best in 28 of 31 evaluations. Inter-cell attention ablation: it only beats the no-inter-cell-attention variant when the **number of unique cells in the set > 32**; validation loss is best at set size 256, reconstruction metrics best at 128. **Genetic perturbation**: on Xaira (X-Atlas:Orion) genetic-perturbation classification, "absolute performance is low for all methods" (the authors attribute this to measurement noise and expression similarity across perturbations). For weak cytokines such as TNF-α and IL-6 across datasets, "pseudo-bulk and DE scores of all methods are close to 0". "Closest/same-cell-type prompt cells" is a very strong baseline; in setting 3, the oracle DonorMean has the best pseudobulk correlation.
- **Code/weights**: github.com/ArcInstitute/stack (code CC BY-NC-SA 4.0); HF `arcinstitute/Stack-Large` (217M parameters, Arc non-commercial model licence).
- **v2 content**: the v2 abstract adds DiseasePert-3M (T cells from 40 donors × 14 diseases, 11 cytokines) and Perturb Sapiens with 892 perturbations (including genetic). Whether the v2 main text contains cross-cell-line CRISPRi ICL results: **(unverified)** (bioRxiv blocked).
- **Relation to our setting**: This is the most direct precedent for "concatenate prompt cells + target controls and use inter-cell attention", and it is Arc's own work. But the evidence comes mainly from strong-effect cytokines/drugs; on weak genetic perturbations, all methods perform poorly. Its module tokenisation (100 tokens per cell) + axial attention is an off-the-shelf way to keep compute under control.

#### A2. Retrieval-Augmented Generation for Predicting Cellular Responses to Gene Perturbation (PT-RAG)
- **Authors**: Andrea G. Di Francesco (Sapienza / ISTI-CNR), Andrea Rubbi (Cambridge / Sanger), Pietro Liò (Cambridge)
- **Date/links**: 2026-03-07, arXiv 2603.07233; Gen² @ ICLR 2026 workshop
- **Method**: Built on STATE. Two stages: (1) retrieve the top-K=32 candidate perturbations by cosine similarity of GenePT embeddings; (2) score triples [h_ctrl; h_pert; h_cxt_k] with a scoring MLP, use Straight-Through Gumbel-Softmax for differentiable sparse selection, add a sparsity regulariser (λ=0.1), and feed the observed responses of the selected perturbations to the generator as context, trained end-to-end jointly with generation.
- **Setting**: Replogle-Nadig (K562/Jurkat/RPE1/HepG2, 2009 single-gene CRISPRi perturbations). Leave-one-cell-type-out; **30% of the target type's perturbations are visible in training** (few-shot, consistent with the STATE protocol), 70% used for validation/testing.
- **Results**: Pearson-DEG: STATE 0.624 / vanilla RAG 0.396 / PT-RAG 0.633; MAE 0.298/0.429/0.295; W2 646.1/1189.5/633.7; Energy 9.41/14.18/9.33. Removing the sparsity regulariser drops Pearson to 0.134; with the same retrieval pool, vanilla RAG (K=32) Pearson 0.351 vs PT-RAG 0.604 (ablation setting). For the same query gene, overlap of selected perturbations across cell types is only ~19%. FLOPs are ~1.7×.
- **Code**: github.com/difra100/PT-RAG_ICLR (paper CC BY 4.0)
- **Relation to our setting**: The closest work to "using related perturbations as context". **Key lesson: injecting context indiscriminately makes things significantly worse** (vanilla RAG is much worse than STATE); sparse, conditional, learnable selection is required. The gain over STATE is small (+0.009 Pearson-DEG). The target line is few-shot, not zero-shot.

#### A3. PerturbMap: Cross-Context Transfer of Single-Cell Perturbation Responses
- **Authors**: Panpan Cui, Yiqi Liu, Wenhao Sun (UCAS / ICT-CAS / HKUST)
- **Date/links**: 2026-07-30, arXiv 2607.28090v1
- **Method**: Builds shared rank-K response coordinates (QR) from the training set only; a recipient-local base predictor (inputs: control-anchor gene correlations, control mean/variance/detection rate); a ridge mapping per source→recipient route (penalty scaled by source energy); calibrates each route's interpolation strength α and reliability ρ on disjoint validation anchors; fuses only routes with ρ>0; optional low-capacity logistic "harm gate".
- **Setting**: Frangieh Perturb-CITE-seq (3 contexts, 200 shared perturbations, 5-fold by perturbation identity), and Jiang et al.'s multi-source Perturb-seq (6 recipients). Predicts conditional means only.
- **Results**: MSE 1.581e-3 vs LowRank base 1.649e-3 (relative −4.1%, 161/200 identities improved); **raw copy (copying the source effect directly) gives only −0.86e-5, with a harm rate of 49.5%**; PerturbMap harm rate 19.5%; only 2.82e-6 away from centralised TokPool.
- **Code**: not stated.
- **Relation to our setting**: Provides quantitative evidence that "borrowing effects across contexts requires reliability weighting and harm control". Note that it needs paired, measured anchor perturbations in the recipient — **not zero-shot**.

#### A4. Tabular Foundation Models Are Competitive Cellular Perturbation Predictors Across Biological Scales
- **Authors**: Palla G., Hillsley A., Kim Y.-J., Royer L.A. (Biohub SF)
- **Date/links**: v3 2026-07-19, bioRxiv 10.64898/2026.06.28.735106 (CC-BY); code github.com/royerlab/tfm-perturbation (MIT)
- **Method**: Uses TabICL/TabPFN as general-purpose in-context regressors. Rows = cells, features = expression (optionally PCA-50), target = perturbation-induced change; labelled source cell types are matched to unlabelled target cell types via nearest neighbours or OT to build the in-context sample set (per the README).
- **Setting**: (1) cell-level ICL across cell types (OpenProblems PBMC, T cells → NK/B/myeloid); (2) pseudobulk prediction on Perturb-seq from 5 cell lines; (3) genome-wide CRISPR in primary CD4+ T cells; (4) zebrafish embryo composition. Compared against PRESAGE, scGPT, scLAMBDA, STACK, Prophet.
- **Results (abstract)**: on par or better across cell types; on pseudobulk "consistently outperforms specialised baselines". Specific numbers **(unverified)** (main text blocked).
- **Relation to our setting**: Shows that "in-context regression" is useful in itself but does not need a biology-specific architecture; also suggests that a TabPFN-style low-dimensional ICL on pseudobulk/PCA should be the first strong baseline for any ICL module.

#### A5. X-Cell: Scaling Causal Perturbation Prediction Across Diverse Cellular Contexts via Diffusion Language Models
- **Authors**: Wang C., Karimzadeh M., …, Wang B. (Xaira Therapeutics)
- **Date/links**: 2026-03-20, bioRxiv 10.64898/2026.03.18.712807 (CC-BY-NC); PDF: cdn.xaira.com/papers/X_CELL_V1_0316_final.pdf
- **Data**: X-Atlas/Pisces, **25.6M CRISPRi genome-wide Perturb-seq cells across 16 contexts** (common cell lines, iPSCs, resting/activated Jurkat, multi-lineage differentiated iPSCs).
- **Method**: Masked discrete diffusion over each cell's gene tokens (part of a control cell's gene values are replaced with perturbed values, and the model learns to reconstruct the full perturbed profile); a cross-attention layer every 3 layers attends to 6 perturbation-gene priors (ESM-2, STRING, GenePT, DepMap, JUMP Cell Painting, gene-identity embeddings); trained on groups of 64 cells (set-level distribution-matching loss). 55M (scGPT initialisation) → X-Cell-Ultra 4.87B. **Test-time adaptation (TTA)**: in the target domain, only NTC cells are used for "control→control" self-supervision, fine-tuning self-attention with an MMD loss only, freezing and skipping cross-attention, with ~200–1024 sets of 64 cells, LR 1e-6–1e-5.
- **Results**: Pearson Δ "up to 5×" over prior SOTA (abstract); STATE scores 0.10 on their benchmark (main-text excerpt). Zero-shot: on iPSC-derived melanocyte precursors (1341 genetic perturbations) and primary CD4+ T cells (2 donors × resting/8h/48h), X-Cell-Ultra beats control mean, perturbation mean, STATE, and X-Cell on Pearson Δ / DE Direction Match / MAE. Scaling: training-loss exponent α=0.32, but the validation-loss exponent is only 0.03; DE Pearson goes from 0.291 (83M) → 0.315 (1.6B).
- **Code/weights**: no release statement for code or weights found in the PDF (release status **(unverified)**).
- **Relation to our setting**: (1) the largest public CRISPRi multi-context dataset (downloadability to be verified); (2) **TTA on target-line controls** is complementary to computing the line embedding from controls only, and can be added at zero cost; (3) its cross-context ability comes from scale + prior cross-attention, not cross-cell ICL.

#### A6. CellMSA: Context Modeling for Single-Cell Representation Learning
- **Authors**: Suyuan Zhao, Minghao Liu, Yizhen Luo, Zaiqing Nie (PharMolix / Tsinghua AIR; affiliations **(unverified)**)
- **Date/links**: 2026-09-30, arXiv 2609.38908; code github.com/PharMolix/CellMSA
- **Method**: For each target cell, retrieves 16 same-batch same-type + 16 cross-batch same-type + 8 similar-type cells to form an "MSA", with a relation-type embedding per row; the CellMSA-Module alternates outer-product mean (aggregating co-variation of gene pairs across rows, complexity O(S·G²·d_p²)) and pair-weighted averaging; the resulting gene-pair representation is injected as an attention bias into the target-cell encoder (GenePairformer).
- **Perturbation evaluation**: filtered Replogle (4 cell lines, 100 perturbations), **unseen perturbations** (not unseen cell types), with STATE-ST downstream. Pearson Δ 0.433 vs STATE-SE 0.353 (the paper reports +8.8% relative to the strongest baseline; the two framings differ, as in the original), PRAUC 0.334, DE Overlap 0.215. (paraphrased from a summarisation tool)
- **Relation to our setting**: The first direct implementation of AlphaFold-style row-column / pair representations for single cells, usable as an "MSA-like" design reference; the G² cost shows that pair representations must be built over a module axis rather than 4096 genes.

#### A7. Predicting cellular responses to perturbation across diverse contexts with STATE
- **Authors**: Adduri A.K., Gautam D., …, Goodarzi H., Roohani Y.H. (Arc Institute)
- **Date/links**: bioRxiv 10.1101/2025.06.26.661135 v2 2025-07-10; published in Cell 2026 (S0092-8674(26)00921-9, main text not read); github.com/ArcInstitute/state
- **Method**: State Embedding (SE, trained on 167M observational cells) + State Transition (ST, a transformer over cell sets, MMD loss, ~100M perturbed cells).
- **Setting**: Replogle-Nadig leave-one-cell-line-out, **with 30% of the target line's perturbations added to training (few-shot)**; there is also a zero-shot setting where "the target line has controls only", used to identify strong perturbations.
- **Results (abstract)**: >30% improvement in discrimination on large datasets; can "identify strong perturbations" in new contexts with no perturbation observations at all. A search summary mentions a >17% gain from SE and Spearman >0.5 in Replogle-Nadig/Parse zero-shot — **(unverified) against the original**.
- **Relation to our setting**: Mainstream baseline; its "zero-shot" capability is limited to identifying strong perturbations, which is precisely not our weak-effect regime.

#### A8. Causal Cellular Context Transfer Learning (C3TL)
- **Authors**: Michael Scholkemper, Sach Mukherjee (DZNE Bonn / Bonn / Cambridge)
- **Date/links**: 2026-03-13, arXiv 2603.13051 (CC-BY)
- **Method**: Models pseudobulk Δ. Perturbation encoding Ê_z(p) = mean_γ Θ_z(x_γp, ψ_γ) over all contexts γ in which p was observed (enforcing context invariance); context encoding Ê_ψ(c) = mean_π Θ_ψ(x_cπ, z_π) over perturbations π observed in context c; the two are summed and decoded. Essentially DeepSets-style cross-context set pooling.
- **Setting**: Replogle (4 folds), Parse (24 contexts), Tahoe (48 cell lines); in the held-out context, **8% train + 2% validation** (few-shot), 90% test.
- **Results** (Table 1, Pearson / MSE): Replogle C3TL 0.491 / 0.016 vs State 0.474 / 0.018 vs Mean 0.370; Parse 0.670 vs 0.608; Tahoe 0.777 vs 0.778. ~30× speed-up, memory 2.1GB vs State 19.8GB.
- **Relation to our setting**: **"Mean pooling across contexts" already matches State**, suggesting that gains from cross-cell/cross-line attention may come mostly from aggregation itself; an ICL module should be compared against "simple pooling". However, its context encoding needs some of the target line's perturbations and is unavailable in the zero-shot case.

#### A9. scPILOT: Predicting Single-Cell Perturbation Responses Across Biological Contexts With a Deep Generative Model Integrating Optimal Transport
- **Authors**: Wang J., Liu Z., Zhang Z., et al. (yongzhuangliulab)
- **Date/links**: Advanced Science, online 2026-08-29; PMC13525497; github.com/yongzhuangliulab/scPILOT (CC BY 4.0)
- **Method**: VAE latent space (discriminator-assisted); first estimates cell-level responses from observed contexts via latent OT, then localises with Leiden and performs query→reference OT transfer, and finally applies adaptive weighting across multiple observed contexts.
- **Results**: cross-cell-type R² 0.945 / MMD² 0.137; cross-patient 0.598 / 0.025; cross-species 0.853 / 0.287; **cell-line benchmark (IFNGR2 knockdown in 6 cancer cell lines) ≈ identity mapping**; the authors describe the perturbation effect as "modest" relative to between-line differences (paraphrased from a summarisation tool).
- **Relation to our setting**: With weak genetic effects + large between-line differences, sophisticated transfer does no better than the identity map — direct counter-evidence.

#### A10. TxPert: Leveraging Biochemical Relationships for OOD Transcriptomic Perturbation Prediction
- **Authors**: Wenkel F., Tu W., Masschelein C., Shirzad H., et al. (Valence Labs / Recursion / UBC)
- **Date/links**: arXiv 2505.14919 (2025-05); Nature Biotechnology 2026 (s41587-026-03113-4, main text blocked)
- **Method**: basal state encoder + GNN perturbation encoding (Exphormer, etc.) over multiple KGs (STRINGdb, GO, PxMap, TxMap, etc.).
- **Cross-line setting**: leave-one-out over K562/RPE1/HepG2/Jurkat, **the test line contributes only controls to training and all its perturbations are unseen (true zero-shot)**; this task **does not use the basal encoder** (the authors say there are too few cell lines to learn a generalisable basal representation). Beats a general baseline and scLAMBDA; numbers are only shown in figures **(unverified)**. Uses split-half "experimental reproducibility" as an upper bound.
- **Code**: GitHub, CC BY-NC-SA 4.0; PxMap/TxMap are proprietary data.
- **Relation to our setting**: The published benchmark most consistent with our core setting (a new line with controls only); its conclusion is that with only 4 cell lines one cannot learn an encoder that "infers context from controls" → supports training a line encoder on more contexts (X-Atlas/Pisces, an H1 hESC CRISPRi Perturb-seq dataset (Arc Institute, 2025), etc.).

#### A11. OCOO-T: A Simple and Scalable Virtual Cell Model
- **Authors**: Danning Jiang, Zheming An, Yalong Zhao, Lipeng Lai (Infevo AI)
- **Date/links**: 2026-06-11, arXiv 2606.12838
- **Method**: Minimalist rectified flow (v-prediction), plain transformer, adaLN conditioning, ESM2 gene-identity embeddings, patching to scale to 18,080 genes, 32 learnable in-context tokens injected at layer 4.
- **Results (paraphrased from a summarisation tool)**: Tahoe PDCorr 0.952 vs STATE 0.535; Replogle 0.437 (on par with STATE), better DE. **Ablation: conditioning on the mean profile of control cells (set size 1–64) performs close to a learnable cell-line embedding** (the latter is stronger on DE discrimination). **MMD regularisation hurts beyond >2000 genes**. The Replogle split appears to follow STATE's few-shot protocol (the training set includes 4 cell lines) — **whether zero-shot or few-shot is (unverified)**.
- **Relation to our setting**: Directly supports the feasibility of computing the line embedding from the control set only; a caution about using MMD on high-dimensional gene axes.

#### A12. STRAND: Sequence-Conditioned Transport for Single-Cell Perturbations
- **Authors**: Boyang Fu, …, Marinka Zitnik (Harvard / Merck / Broad, etc.)
- **Date/links**: 2026-02-10, arXiv 2602.10156 (CC-BY)
- **Method**: Represents perturbations by the regulatory DNA sequence at the target site (DNA foundation model), mapped into an RNA latent, conditioning a control→perturbed latent bridge matching (OT pairing).
- **Setting**: PerturbQA splits; K562/Jurkat/RPE1 (HepG2 excluded because the median number of cells per perturbation is only 45 and DE statistical power is low); includes cross-line "Combined" data and transfer to unseen cell lines.
- **Results**: in the low-sample regime, discrimination score 33% higher than STATE (0.72 vs 0.54); on unseen cell lines, ΔPearson up to +0.14 (abstract). The authors **deliberately exclude metrics that use Wilcoxon DE as hard labels** (citing Squair 2021), arguing they have severe false positives.
- **Relation to our setting**: CRISPRi is TSS/locus-specific; representing perturbations by sequence can explain differences between different sgRNAs for the same gene. When the same sgRNA set is used across lines, locus information helps cross-line consistency.

#### A13. pertTF: context-aware AI modeling for genome-scale and cross-system perturbation prediction
- **Authors**: Su Y., Liu D., Menon V., …, Huangfu D., Li W. (UMaryland / MSKCC, etc.)
- **Date/links**: 2026-03-16, bioRxiv 10.64898/2026.03.12.711379; PMC13015719; github.com/davidliwei/pertTF (MIT), HF weililab
- **Method/setting**: A transformer trained on 30 KOs × 14 pancreatic differentiation cell types; the "new context" setting makes WT cells of a held-out type visible during training but not its perturbed cells (isomorphic to our core setting), and there is also a joint unseen-gene × unseen-type setting. Numbers **not individually verified**.
- **Relation to our setting**: Another example of an isomorphic setting, but small-scale and with in-house data.

#### A14. AlphaCell (Tongji) / Lingshu-Cell (Alibaba DAMO)
- AlphaCell: 2026-03-05, bioRxiv 10.64898/2026.03.02.709176; whole-transcriptome latent + OT-CFM; claims zero-shot prediction in "completely unseen cellular contexts". **Abstract only; results (unverified)**.
- Lingshu-Cell: 2026-03-26, arXiv 2603.25240 (ICML 2026); masked discrete diffusion over ~18k genes; joint embedding of cell identity × perturbation; claims leading results on an H1 hESC CRISPRi Perturb-seq dataset (Arc Institute, 2025). **Abstract only**.
- **Relation to our setting**: Low-priority references.

### 2.B Flow matching / OT / diffusion generators

#### B1. PRiMeFlow (Altos Labs)
- **Authors**: Zichao Yan, Yan Wu, Mica Xu Ji, …, Rory Stark
- **Date/links**: 2026-04-15, arXiv 2604.13986; github.com/altoslabs/primeflow
- **Method**: CFM directly in full-gene log1p space, with a U-Net velocity field (based on guided-diffusion); **the source distribution is a standard Gaussian (not control cells)**, independent coupling, conditioning = sum of perturbations + one-hot covariates, CFG (1–5) at inference.
- **Results**: on PerturBench covariate transfer (Srivatsan20), MMD 0.13 and DEG recall 0.26 are both best; FM-PCA's DEG recall is only 0.042; replacing the U-Net with an MLP degrades badly; OT pairing in PCA gives no clear benefit. **Worse than CPA/SAMS-VAE on pseudobulk metrics (LogFC, RMSE)**. (paraphrased from a summarisation tool)
- **Relation to our setting**: Supports end-to-end FM in gene space; but its advantage lies in higher-order distributional moments, whereas many common evaluation metrics for perturbation prediction (log-fold-change error, direction agreement, DE-gene-set overlap such as Jaccard) are pseudobulk/DE-oriented → a control-anchored design is more compatible with pseudobulk accuracy.

#### B2. scDFM: Distributional Flow Matching for Robust Single-Cell Perturbation Prediction (ICLR 2026)
- **Authors**: Chenglei Yu, Chuanrui Wang, Bangyan Liao, Tailin Wu (Westlake University / Zhejiang University)
- **Date/links**: 2026-02-06, arXiv 2602.07103; github.com/AI4Science-WestlakeU/scDFM (CC BY 4.0)
- **Method**: CFM + multi-kernel MMD regularisation; PAD-Transformer (co-expression kNN graph with k=30 as attention mask; differential attention, applied twice to the control representation as self- and cross-attention).
- **Results**: Norman additive MSE 0.00315 vs CellFlow 0.00392 (−19.6%); Norman holdout discrimination 0.9189. Tested only within K562/A549, **no unseen cell types**. Removing MMD causes the largest drop.
- **Relation to our setting**: MMD helps in their setting (~5000 genes), contradicting OCOO-T (MMD hurts beyond 2000 genes); this needs our own ablation on a 4096-gene axis.

#### B3. GeneGeoFlow: Control-Anchored Residual Flow Matching Conditioned on Gene Geometry
- **Authors**: Li, Chi, Song, Zhang, et al.
- **Date/links**: 2026-08-07, arXiv 2608.06824
- **Method**: Learns a residual velocity field starting from (noised) control cells; uses multi-scale spectral coordinates from GO + control co-expression graphs as gene geometry, selected by perturbation-conditioned Spectral Scale Router / Graph Source Router; within-condition OT pairing; a **Δ-correlation objective** aligns the direction of predicted and true condition-level shifts.
- **Results**: Norman strict 5-fold Pearson Δ 0.8153 (graph-free 0.1450; this number looks suspicious, as in the original); ComboSciPlex 0.9088 vs scDFM 0.8933; routing +0.0357; explicit graph propagation gives no consistent gain; shuffled coordinates give 0.4387. **Unseen cell types not tested**. Code not stated.
- **Relation to our setting**: Isomorphic to the control-anchored transport approach we explored earlier; the "Δ-direction correlation loss" is directly relevant to direction-type evaluation metrics (direction fidelity, cosine-based perturbation discrimination).

#### B4. scBIG: Beyond Independent Genes — Module-Inductive Representations
- **Authors**: Jiafa Ruan, Ruijie Quan, Liyang Xu, Zongxin Yang, Yi Yang (Zhejiang University, etc.)
- **Date/links**: v2 2026-06-16, arXiv 2602.04901; github.com/ttruan2426-dot/scBIG
- **Method**: Gene-Relation Clustering (fusing foundation-model semantics + high-confidence relations, inducing gene modules via OT) → Gene-Cluster-Aware Encoder (inducing-point bottleneck attention modelling inter-module interactions) → latent CFM; Cluster Correlation Alignment + Pathway-informed OT regularisation.
- **Results**: +6.7% on average over the strongest baseline; tests unseen/combinatorial perturbations on Norman, K562, RPE1 (**each cell line trained separately, not cross-line**).
- **Relation to our setting**: Supports a "co-response module axis"; modules can serve as the unit of ICL tokens.

#### B5. PerturbCellRL: Verifier-Guided RL for Single-Cell Perturbation Prediction
- **Authors**: Dongxia Wu, Mingyu Li, …, Emma Lundberg, Serena Yeung-Levy, Emily B. Fox (Stanford / PKU)
- **Date/links**: 2026-06-26, arXiv 2606.27752
- **Method**: Uses 4 cell-level verifiers (Pearson top-k, RMSE top-k, DE Spearman, pathway activity) as rewards for RL post-training of a pre-trained FM generator (scDFM-like); also usable for test-time scaling. Datasets include Norman, ComboSciPlex, and an H1 hESC CRISPRi Perturb-seq dataset (Arc Institute, 2025).
- **Relation to our setting**: An option for "metric-targeted post-training" on top of a control-anchored flow-matching model; orthogonal to ICL.

#### B6. CellFlow (2025-04) and CFM-GP (2025-08) — necessary 2025 background
- CellFlow: Klein D., …, Theis F.J. (Helmholtz Munich), bioRxiv 10.1101/2025.04.11.648220 (CC-BY-NC-ND); FM for multimodal conditions including cytokines, drugs, KOs, organoid protocols, etc.
- CFM-GP: Abir, Dip, Zhang (Virginia Tech), arXiv 2508.08312; cell-type-conditioned CFM, cross-cell-type in 5 cytokine/infection/drug scenarios (not cross-line CRISPR).
- **Relation to our setting**: The precedents for cross-context FM are all strong-effect stimuli/drugs.

### 2.C Knowledge graph / pathway / module priors

#### C1. PRESAGE (Genentech, 2025-06)
- Littman R., …, Regev A., Hütter J.-C.; bioRxiv 10.1101/2025.06.03.657653 (CC-BY-NC-ND)
- Aggregates gene embeddings from 40 knowledge sources (node2vec for graphs, PCA for tables); an attention model predicts unseen single-gene perturbations. Ablation conclusions: **choice of knowledge sources matters more than architectural complexity; "cross-system Perturb-seq data" is a particularly predictive gene embedding**; performance saturates quickly with training-set size, and the authors recommend "sparse measurement across many systems".
- **Relation to our setting**: The most direct implication is to treat "same-target Perturb-seq effects in other cell lines" as a gene prior (a static embedding), without necessarily using ICL.

#### C2. Knowledge Graphs and Reasoning LLMs for Finding Simple Yet Effective Transcriptomic Perturbation Predictors
- Jake Fawkes, Liam Hodgson, Jason Hartford (Valence Labs, etc.; **affiliations (unverified)**); 2026-06-07, arXiv 2606.08816
- **KG-based kNN (the simplest model) is already "highly competitive"**; RL-optimised reasoning LLMs can further match SOTA and improve OOD perturbations. No code seen.
- **Relation to our setting**: When selecting "pathway-related perturbations" for ICL, KG-kNN is the natural retriever and baseline.

#### C3. AdaPert: Learning Adaptive Perturbation-Conditioned Contexts (ICML 2026)
- Piao Y., …, Park C., Ahn S.; arXiv 2602.18885 (v2 2026-07-05)
- Uses Gumbel-Softmax to select sparse, perturbation-specific subgraphs from a KG; a Huber robust loss on non-DEGs + an alignment loss on DEGs to counter mean-collapse. K562: Pearson-Δ 0.619 vs TxPert 0.580; DES@50 0.263 vs 0.220; **small-effect DEG recall 11.5% vs 9.0%**. Includes a cross-line study (4 cell lines, numbers in the appendix, **(unverified)**). No code seen.
- **Relation to our setting**: Loss design for weak-effect perturbations (suppressing spurious changes in non-responsive genes) is directly relevant.

#### C4. Stable-Shift (BCB '26) / PerturbGraph (bioRxiv 2026-03) / CisTransCell (2026-06)
- Stable-Shift (Dip & Zhang, arXiv 2606.24940): a low-rank response basis from training perturbations + a GCN over STRING/GO/control statistics predicting coordinates for unseen genes; K562 cosine 0.592 vs GEARS 0.569.
- PerturbGraph (Dip & Zhang, bioRxiv 10.64898/2026.03.23.713780, github.com/Sajib-006/PerturbGraph): latent program space + GNN; cosine +6% over tree models.
- CisTransCell (Zhang W., et al., Leibniz Hannover, arXiv 2606.13713): regulatory-sequence + coding-sequence priors, FiLM conditioning; compared with PRESCRIBE, tests unseen perturbations only.
- **Relation to our setting**: All are limited to unseen perturbations within a single cell line; "low-rank response basis + prior-predicted coordinates" can serve as a zero-shot baseline.

### 2.D Benchmarks / critical evaluations / metrics

#### D1. Virtual Cells Need Context, Not Just Scale (position paper)
- Dibaeinia P., Babu S., Knudson M., ElSheikh A., Wen Y., Liu H., Perera J., Khan A.A. (Biohub / UChicago / Northwestern); 2026-02-09, bioRxiv 10.64898/2026.02.04.703804 (CC-BY); PMC12919078 (full text read)
- Argument: context is an effect modifier, f_c changes with c; under causal transportability, it is not identifiable without target-domain interventions or invariance assumptions. Empirics: Zhu 2025 primary CD4+ T cells (22M cells, ~12k knockdowns, 4 donors × 3 time points = 12 contexts), cross-context transfer with scLDM. Correlation-Δ / discrimination (median >0.9) look good, but **DEG precision ~0.67 and recall ~0.09**; DEG-F1 correlates only weakly with aggregate metrics (r=0.26, −0.05). **As context coverage increases from ≤3 to 8, DEG-F1 rises from <0.1 to ~0.19**; correlation with cell count is only r=0.11. At comparable cell counts, the high-diversity group is significantly better (p=0.002).
- **Relation to our setting**: Training should maximise cell-line / context coverage rather than cell count; in the weak-effect regime, DE recall is the core difficulty.

#### D2. Mechanisms Matter: Transportability of Cellular Perturbation Effects
- Qi S.-a., Chapfuwa P. (Microsoft); v2 2026-07-02, bioRxiv 10.64898/2026.05.08.723625 (CC-BY)
- A semi-synthetic causal simulator (tunable mechanism differences) + Vendi diversity score to diagnose mode collapse; 4 DL models + 6 simple baselines: **under cross-context splits, performance often drops to the level of simple baselines; even on synthetic data with fully known causal structure, no model generalises across contexts with different mechanisms**. (abstract only)
- **Relation to our setting**: ICL cannot infer unseen mechanism differences out of nothing; it can only help when the mechanisms of the "related contexts" are close to those of the target line → need a way to judge which source lines are mechanistically close to the target (reliability weighting).

#### D3. Can We Trust In-Distribution Success? Locked Evaluation … CRISPRi Perturbation Prediction
- Mehrdad Shoeibi, Niloofar Yousefi (UCF); arXiv 2608.00152 v3 (2026-08-31)
- Frozen Geneformer representations are informative in-distribution on an H1 hESC CRISPRi Perturb-seq dataset (Arc Institute, 2025) (ΔR² +0.1645 vs random features), but **zero-shot transfer to two Replogle screens is negatively correlated** (ρ=−0.139, −0.267); **effect-size endpoints on that dataset are strongly correlated with cell counts** (counts alone R²=0.4325 vs 4 magnitude scalars 0.2589). Cites Ahlmann-Eltze 2025: linear perturbation embeddings pre-trained on Replogle K562 transfer to RPE1 (and vice versa).
- **Relation to our setting**: Evaluation must control for cells per perturbation (which varies widely across training source screens); endpoints must be harmonised across screens.

#### D4. Score Distributions, Not Cells: Evaluating Single-Cell Perturbations Under Class Overlap
- Marrakchi Y., D'Ascenzo D., Cultrera di Montesano S.; 2026-07-06, arXiv 2607.04595
- Single-cell distributions of weak perturbations overlap heavily with controls; single-cell classification F1 is ~0.2–0.3; averaging classifier probabilities over the whole population and ranking (CDS) gives perturbation-level accuracy of 0.976–1.000. On an H1 hESC CRISPRi Perturb-seq dataset (Arc Institute, 2025; 301 labels), CDS rank-1 is 0.86 vs ~0.70 for the best PDS variant; with only 10% of cells, PDS (especially L1) degrades markedly. (paraphrased from a summarisation tool)
- **Relation to our setting**: Directly shows that under weak effects, the single-cell-level "cross-cell attention" signal is weak and the information lives mainly in population statistics → ICL tokens should be at the population / pseudobulk / module level.

#### D5. Evaluating Single-Cell Perturbation Response Models Is Far from Straightforward
- Heidari M., Karimpour M., Srivatsa S., Montazeri H. (U Tehran); 2026-02-17, bioRxiv 10.64898/2026.02.14.705879
- Correlation metrics, Wasserstein, and Energy are affected by scale, sparsity, and dimensionality; Wasserstein fails under high-dimensional variance scaling, and Energy can miss broken gene-gene dependencies; on chemical perturbation data, complex models often underperform simple baselines. (abstract only)

#### D6. Deep Learning-Based Genetic Perturbation Models Do Outperform Uninformative Baselines on Well-Calibrated Metrics (2025-10, necessary background)
- Miller H.E., Mejia G.M., … (Shift Bioscience); bioRxiv 10.1101/2025.10.20.683304 (CC-BY)
- Uses interpolated duplicates as positive controls and the dynamic range fraction as a calibration measure; 14 datasets × 13 metrics: **MSE and Pearson(Δctrl) are poorly calibrated; weighted or rank-based metrics are well calibrated**; under well-calibrated metrics, DL outperforms mean / control / linear baselines.

#### D7. Foundation Models Improve Perturbation Response Prediction (GenBio AI)
- Cole E., …, Bar-Joseph Z., Xing E.P.; 2026-02-19, bioRxiv 10.64898/2026.02.18.706454; github.com/genbio-ai/foundation-models-perturbation (GenBio community licence)
- 600+ models: some FMs significantly beat simple baselines, multi-FM attention fusion is better still, and with enough data performance approaches the ceiling. The README says latent diffusion / flow matching and similar methods lag behind simple methods built on FM embeddings, and that chemical perturbations are harder than KOs — specific numbers for cross-cell-line tasks **(unverified)**.

#### D8. VCBench: Benchmarking virtual cell models for in-the-wild perturbation response
- Xinjie Mao, et al. (Shanghai AI Lab); 2026-04-30, arXiv 2604.27646; github.com/maoxinjie/VCBench
- Embedding-based unseen-cell splits are much harder than random splits; **on unseen-cell tasks, linear additive / BioLORD / scLAMBDA are the most stable**; for unseen perturbations, linear methods are insufficient; naively merging datasets can hurt; rankings vary widely across metrics. (paraphrased from a summarisation tool)

#### D9. What Makes a Representation Good for Single-Cell Perturbation Prediction?
- Jiang, et al.; 2026-05, arXiv 2605.19343
- FM representations are less linearly decodable for perturbation labels than PCA; advocates contrastive alignment that separates the "perturbation-invariant background" from the "sparse perturbation effect". (paraphrased from a summarisation tool; low priority)

#### D10. TRADE: Transcriptome-wide analysis of differential expression in perturbation atlases (Nadig et al., Nat Genet 2025-04, necessary background)
- Estimates "transcriptome-wide impact": **a typical gene perturbation affects ~45 genes, a typical essential gene >500**; many effects are undetectable with standard DE; **cross-cell-type consistency of perturbation effects is "moderate"**. (Europe PMC abstract)
- **Relation to our setting**: When the targets of interest are dominated by non-essential genes, there are very few true DEGs per target; the ceiling for cross-line borrowing is bounded by "moderate consistency".

---

## 3. Synthesis

### (a) Established findings

1. **Cross-context generalisation, not model capacity, is the main bottleneck.** Under cross-context splits, all kinds of models often degrade to simple baselines (Mechanisms Matter, VCBench, Virtual Cells Need Context); even on synthetic data with known causal structure, there is no generalisation across contexts with different mechanisms.
2. **Context diversity takes priority over cell count**: on T-cell data, increasing training context coverage from ≤3 to 8 roughly doubles DEG-F1, essentially independently of cell count (VCNC); PRESAGE also recommends sparse measurement across many systems.
3. **Weak effects are the norm**: non-essential gene perturbations affect only ~45 genes on average, with moderate cross-cell-type consistency (TRADE); DEG recall under cross-context transfer is extremely low (~0.09, VCNC).
4. **Set / population-level aggregation works**: STATE (sets + MMD), Stack (inter-cell attention, beneficial only beyond >32 unique cells), C3TL (cross-context mean pooling matches State), CDS (discrimination after population averaging is far stronger than single-cell discrimination).
5. **Characterising a cell line from its control set alone is feasible**: OCOO-T conditions on the mean control profile and performs close to a learnable cell-line embedding; X-Cell's TTA on target-line NTCs improves zero-shot performance. Conversely, TxPert notes that 4 cell lines are not enough to learn a generalisable basal encoder.
6. **FM generator design**: FM in gene space + a strong backbone (U-Net / transformer) helps distributional fit and DEG recall (PRiMeFlow), but is not superior on pseudobulk LFC; a control-anchored residual + Δ-direction loss (GeneGeoFlow) fits direction-type evaluation metrics better; the effect of MMD depends on dimensionality (positive in scDFM, negative beyond >2000 genes in OCOO-T).
7. **Priors > architecture**: KG-kNN is already strong (Fawkes et al.); PRESAGE finds "cross-system Perturb-seq effects" to be the most useful gene embedding.
8. **Metric pitfalls**: MSE and Pearson(Δctrl) are poorly calibrated; L1-PDS degrades markedly with fewer cells (Score Distributions); effect-size endpoints are entangled with cell counts; Wilcoxon DE as hard labels yields many false positives (STRAND excludes such metrics for this reason), so DE-based metrics such as Wilcoxon significance overlap or DE-gene-set Jaccard should be interpreted with care in the weak-effect regime.

### (b) Remaining gaps: "zero-shot, genetic, new line ← pathway-related perturbed cells from other lines as context"

- **No paper was found that directly evaluates this setting.** The closest works each miss one piece:
  - PT-RAG: retrieves **related perturbations** as context ✔, genetic ✔, but the target line is few-shot (30% visible) ✘, and it operates at STATE's pseudobulk / set level, not cell-level cross-attention.
  - Stack: cell-level prompt/query ICL ✔, zero-shot new context ✔, but the prompt consists of samples of the **same condition** in other cell types; no positive results for genetic CRISPRi across lines (all methods are weak on Xaira genetic perturbations), and whether v2 fills this gap is **(unverified)**.
  - PerturbMap / TabPFN / scPILOT: cross-context transfer of the same perturbation ✔, but PerturbMap needs anchors in the target context; scPILOT ≈ identity on the cell-line genetic task.
  - CellMSA: MSA-style context ✔, but context cells are used for representation learning, and it tests unseen perturbations, not new lines.
  - TxPert / pertTF / X-Cell: zero-shot new context ✔, but relying on priors / scale, not ICL.
- In our setting, the target line has **no perturbed cells at all**, so the ICL context can only be perturbed cells from "{same target gene, pathway neighbours} × {other lines}" + target-line controls. This amounts to merging PerturbMap / TabPFN-style "same-target cross-line transfer" with PT-RAG-style "related-perturbation retrieval" in one cross-cell attention model — **there is no precedent in the literature, and no ablation evidence**.
- **Unresolved identification problem**: with controls only, how can "the target line's mechanism f_c" be inferred? Following the transportability arguments in VCNC / Mechanisms Matter, this requires invariance assumptions; the only available proxies are: control-set embeddings (OCOO-T; TxPert's failure), TTA on controls (X-Cell), and "weighting source lines by control similarity to the target line" (no published evaluation found).

### (c) Is cross-cell attention useful for weak-effect / noisy genetic perturbations? The evidence

**For**
- Stack: the inter-cell attention ablation shows gains on the validation set, but only with >32 unique cells; the authors conclude that "aggregating context is crucial for subtle states".
- STATE: set modelling yields >30% improvement in discrimination on large datasets.
- PT-RAG: under cross-cell-type few-shot, learnable sparse related-perturbation context gives a small gain over STATE (Pearson-DEG 0.624→0.633, W2 646→634, FDR<0.01).
- CellMSA: MSA-style context gives +8.8% Pearson Δ (unseen perturbations).
- PerturbMap: reliability-weighted cross-context borrowing reduces MSE by 4.1% and lowers the harm rate from 49.5% (raw copy) to 19.5%.

**Against / caution**
- PT-RAG: **indiscriminate context is clearly harmful** (vanilla RAG Pearson-DEG 0.396 vs STATE 0.624; 0.134 after removing the sparsity constraint).
- PerturbMap: copying source effects directly gives almost no gain and harms about half of the perturbations; the gap to simple centralised TokPool is tiny (2.8e-6).
- scPILOT: ≈ identity on the cell-line IFNGR2 knockdown task (effects too small relative to between-line differences).
- Stack: all methods score low on genetic-perturbation classification; for weak cytokines (TNF-α / IL-6) across datasets all methods are near 0; closest-cell-type prompts and DonorMean are very strong baselines.
- C3TL: simple mean pooling already matches State → the extra benefit of attention may be limited.
- Score Distributions: under weak effects, single cells are nearly indistinguishable and information concentrates in population statistics → cell-level cross-cell attention gets almost no signal per single cell; its main role is implicit averaging.
- Mechanisms Matter / VCNC: no transfer across contexts with different mechanisms, regardless of architecture.

**Net assessment**: for weak-effect genetic perturbations, **there is no evidence that cell-level cross-cell attention beats "population / module-level aggregation + selective weighting"**. The positive evidence comes entirely from strong effects (cytokines, drugs) or few-shot settings; the negative evidence consistently points to "unfiltered context is harmful" and "between-line mechanism differences set the ceiling". Plausible sources of gain are: (1) denoising by averaging same-target effects across multiple lines; (2) learning when to trust which source / neighbour (reliability weighting). Neither necessarily requires cell-level attention.

### (d) Concrete recommendations (for extending a control-anchored flow-matching approach)

1. **Establish a "transfer-baseline ceiling" before deciding on ICL** (lowest cost, 1–2 days): on leave-one-cell-line-out (HepG2 / Jurkat / RPE1 / K562 / X-Atlas lines / the H1 hESC CRISPRi Perturb-seq dataset (Arc Institute, 2025), etc.), measure (i) same-target cross-line mean transfer, (ii) mean transfer + magnitude restoration (norm-restore), (iii) effect averaging over KG-kNN neighbours, (iv) the noise ceiling from split-half. The ICL module is only worth keeping if it clearly exceeds (ii)+(iii) (this is the kill criterion).
2. **If doing ICL, use population / module-level tokens, not raw cells**: each (source line, perturbation) yields one or a few "effect tokens" (pseudobulk Δ projected onto co-response modules of the 4096-gene axis + the line's control embedding + the perturbation prior). Target-line control cells (the flow state of the control-anchored model) read these tokens via cross-attention. Rationale: TRADE / CDS / Stack's genetic results all show the single-cell signal is too weak; Stack's 100 module tokens show that module-level tokens are viable; CellMSA's G² cost shows pair representations cannot be built at gene level.
3. **Make selection and reliability explicit in attention**: PT-RAG-style sparse Gumbel selection (with a sparsity regulariser) + PerturbMap-style route reliability (estimating each source line's trustworthiness on validation anchors); treat same-target tokens and "pathway-neighbour" tokens separately (different type embeddings, like CellMSA's relation types).
4. **Residual + zero-initialised gate, guaranteeing no regression versus the base model**: the ICL branch only adds a residual to the velocity field, `v = v_base + g·Δv_ICL`, with g initialised to 0; compute the per-target harm rate (fraction of targets that get worse relative to the base control-anchored model) as one of the primary acceptance metrics (the PerturbMap lesson).
5. **Train in episodes that simulate testing**: each episode holds out one line, whose perturbed cells never enter the prompt; prompts come only from other lines (same target + KG/module neighbours). Maximise line coverage (VCNC: diversity over quantity), and downsample to align cells per perturbation (the locked-evaluation paper's finding: effect-size endpoints are entangled with cell counts).
6. **Attention implementation**: if insisting on cell-level attention, use Stack-style axial attention (intra-cell module-token attention + inter-cell attention), with ≥32 unique cells in the set (Stack's threshold); use FlexAttention block masks to implement "target-line controls → prompt" one-way reading, with prompts blocked by (line, perturbation), avoiding the full O((N_ctrl+N_prompt)²) cost and information leakage; use varlen packing for varying numbers of sources.
7. **Zero-cost add-ons**: X-Cell-style TTA (target-line NTCs only, control→control, MMD, tuning only the line encoder / self-attention, freezing the perturbation path); OCOO-T's findings support continuing to compute the line embedding from the control set only.
8. **Align losses with the evaluation goals for weak-effect perturbations**: add a Δ-direction correlation loss (GeneGeoFlow) and Huber suppression on non-DEGs (AdaPert, +28% small-effect DEG recall); internally report well-calibrated metrics (rank / weighted), harm rate, and results stratified by cell count; note that distance-based discrimination scores (e.g. L1-PDS) are sensitive to cell count (Score Distributions), and the effects of magnitude restoration on discrimination and log-fold-change error metrics should be examined separately.
9. **A prior-based alternative (to compare against ICL)**: following PRESAGE, feed "same-target Perturb-seq effects from other lines" directly as a static gene embedding (one vector per target gene) into the perturbation conditioning of the base flow-matching model — effectively "amortised ICL". The cost is far lower than cross-cell attention, and it should be the direct control arm for ICL.

---

## Sources

- Stack preprint: https://www.biorxiv.org/content/10.64898/2026.01.09.698608v2 (PMC full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC12803207/ ; Europe PMC API: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12803207/fullTextXML ; bioRxiv API: https://api.biorxiv.org/details/biorxiv/10.64898/2026.01.09.698608 )
- Stack Arc news: https://arcinstitute.org/news/foundation-model-stack
- Stack-Large HF: https://huggingface.co/arcinstitute/Stack-Large
- Stack GitHub: https://github.com/ArcInstitute/stack
- PT-RAG: https://arxiv.org/html/2603.07233 ; https://github.com/difra100/PT-RAG_ICLR
- CellMSA: https://arxiv.org/html/2609.38908
- PerturbMap: https://arxiv.org/html/2607.28090v1
- Tabular FMs: https://api.biorxiv.org/details/biorxiv/10.64898/2026.06.28.735106 ; https://github.com/royerlab/tfm-perturbation
- X-Cell: https://api.biorxiv.org/details/biorxiv/10.64898/2026.03.18.712807 ; https://www.cdn.xaira.com/papers/X_CELL_V1_0316_final.pdf
- STATE: https://api.biorxiv.org/details/biorxiv/10.1101/2025.06.26.661135 ; https://www.cell.com/cell/fulltext/S0092-8674(26)00921-9 ; https://github.com/ArcInstitute/state
- C3TL: https://arxiv.org/pdf/2603.13051
- scPILOT: https://pmc.ncbi.nlm.nih.gov/articles/PMC13525497/ ; https://pubmed.ncbi.nlm.nih.gov/42667121/
- TxPert: https://arxiv.org/html/2505.14919v1 ; https://www.nature.com/articles/s41587-026-03113-4
- OCOO-T: https://arxiv.org/html/2606.12838v1
- STRAND: https://arxiv.org/pdf/2602.10156
- pertTF: https://pmc.ncbi.nlm.nih.gov/articles/PMC13015719/
- AlphaCell: https://api.biorxiv.org/details/biorxiv/10.64898/2026.03.02.709176
- Lingshu-Cell: https://arxiv.org/abs/2603.25240
- PRiMeFlow: https://arxiv.org/html/2604.13986v1 ; https://github.com/altoslabs/primeflow
- scDFM: https://arxiv.org/html/2602.07103v1 ; https://github.com/AI4Science-WestlakeU/scDFM
- GeneGeoFlow: https://arxiv.org/html/2608.06824
- scBIG: https://arxiv.org/pdf/2602.04901 ; https://github.com/ttruan2426-dot/scBIG
- PerturbCellRL: https://arxiv.org/pdf/2606.27752
- CellFlow: https://api.biorxiv.org/details/biorxiv/10.1101/2025.04.11.648220
- CFM-GP: https://arxiv.org/abs/2508.08312
- PRESAGE: https://api.biorxiv.org/details/biorxiv/10.1101/2025.06.03.657653
- KG + reasoning LLMs: https://arxiv.org/abs/2606.08816
- AdaPert: https://arxiv.org/abs/2602.18885 ; https://arxiv.org/html/2602.18885v2
- Stable-Shift: https://arxiv.org/abs/2606.24940
- PerturbGraph: https://api.biorxiv.org/details/biorxiv/10.64898/2026.03.23.713780
- CisTransCell: https://arxiv.org/pdf/2606.13713
- Virtual Cells Need Context: https://www.biorxiv.org/content/10.64898/2026.02.04.703804v1 ; https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12919078/fullTextXML
- Mechanisms Matter: https://api.biorxiv.org/details/biorxiv/10.64898/2026.05.08.723625
- Locked evaluation CRISPRi: https://arxiv.org/pdf/2608.00152
- Score Distributions, Not Cells: https://arxiv.org/html/2607.04595
- Evaluating … Far from Straightforward: https://api.biorxiv.org/details/biorxiv/10.64898/2026.02.14.705879
- Well-calibrated metrics (Shift Bio): https://api.biorxiv.org/details/biorxiv/10.1101/2025.10.20.683304
- Foundation Models Improve (GenBio): https://api.biorxiv.org/details/biorxiv/10.64898/2026.02.18.706454 ; https://github.com/genbio-ai/foundation-models-perturbation
- VCBench: https://arxiv.org/html/2604.27646v1
- What Makes a Representation Good: https://arxiv.org/html/2605.19343v1
- TRADE: https://www.nature.com/articles/s41588-025-02169-3 (Europe PMC abstract PMC13063516)
