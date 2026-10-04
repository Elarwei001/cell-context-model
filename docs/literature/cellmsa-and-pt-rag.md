# CellMSA and PT-RAG: A Close Reading and Implications for Zero-Shot, Cross-Cell-Line In-Context Perturbation Prediction

> **Motivating setting.** Predicting perturbation responses in cellular contexts never seen perturbed: e.g. CRISPRi knockdowns in a cell line for which only unperturbed control cells are available at training time. Many perturbations in this setting have weak effects, and the goal is to predict the distribution of perturbed cells (not just a mean shift) for each target.
>
> **Design under consideration.** Assemble a set of cells (controls from the new cell line + perturbed cells from other cell lines for pathway-related targets) and make in-context predictions with bidirectional cross-cell attention.
>
> **Reading date:** 2026-10-04. All numbers are taken from the arXiv HTML full texts and the official GitHub code. **(unverified)** marks content that is not stated explicitly in the paper or code and can only be inferred. **(our analysis)** marks our own inferences or calculations, not claims made by the paper authors.

---

## 0. One-page summary (TL;DR)

| | CellMSA (arXiv 2609.38908, NeurIPS 2026) | PT-RAG (arXiv 2603.07233, ICLR 2026 Gen² workshop) |
|---|---|---|
| What it is | A single-cell **representation-learning** foundation model (47.13M parameters). "Context cells" are compressed into a gene×gene pair representation, which is then injected as an attention bias into a gene-token Transformer that encodes **the target cell alone** | A "retrieval" module added on top of STATE-ST. What is retrieved is **only the GenePT text embeddings of neighbouring target genes**, with no expression data at all; Gumbel-Softmax decides which of them enter a weighted sum |
| Cross-cell attention? | **No.** Context rows (cells) interact only through an outer-product **mean** (invariant to row order, linear cost in the number of rows); there is no MSA-Transformer-style column attention | The backbone is STATE's cell-set Transformer (sets of 64 cells), but the retrieved content is target embeddings, not cells |
| Perturbation task setting | 4 Replogle cell lines, 100 perturbations; 45% of HepG2 perturbations held out for testing, while **the other 50% of HepG2 perturbations and all perturbations of the other 3 cell lines are used for training**, so this is not a zero-shot cell-line setting | Replogle-Nadig, 4 cell lines, 2,009 perturbations; **30% of perturbations in the target cell line are visible** (few-shot), so this is not a zero-shot cell-line setting either |
| Key numbers | Pearson Δ: CellMSA+STATE-ST 0.433, raw expression+STATE-ST 0.398 (+8.8%), STATE-SE 0.353, Stack 0.358; single fixed split, single run | Pearson DEG: PT-RAG 0.633, STATE+GenePT 0.631 (difference not significant, pFDR=0.590), STATE 0.624; Vanilla RAG 0.396 (clear degradation) |
| Most useful takeaways for this setting | (1) Add a "relation-type embedding" to context rows; (2) **diversity** of context (cross-batch, related cell types) matters more than **quantity**, with saturation at about 40 cells; (3) mask the same genes in context and target simultaneously so the model cannot copy the answer; (4) use an outer-product mean to compress context into a gene-pair prior, at a cost linear in the number of context rows | (1) Representing the target with GenePT already captures most of the gain; (2) **context injection needs a residual or gate**: Vanilla RAG **replaces** the ctrl+pert representation with the cross-attention output, which is likely the main reason for its degradation (our analysis); (3) without a sparsity constraint, all 32 candidates get selected and performance collapses to Pearson 0.134 |
| Caveats | The perturbation experiment uses a single split and only the 100 perturbations with the most cells (biased toward strong effects); the gene axis is truncated to the target cell's **non-zero genes**, so a target gene knocked down to 0 disappears from the input | The reported evidence for "cell-type-aware retrieval" (cross-cell-line Jaccard ≈0.19) matches the 0.191 expected under random selection from a shared 32-candidate pool (our analysis, see §B.7); MMD is identical across all models, which may indicate an evaluation issue; some numbers in the main text differ from those in the tables |

**In one sentence:** neither paper has been validated on "zero-shot new cell line + weak effects". CellMSA offers an engineering template for how to organise context, how to label the role of each row, and how to summarise cross-cell statistics cheaply. PT-RAG is most informative through its negative results (context injection without a residual degrades; if "retrieval" only retrieves embeddings, its value is roughly that of a better target embedding).

---

## Background: what a protein MSA is, how MSA Transformer's row/column (axial) attention works, and how CellMSA maps this to cells × genes

### What is an MSA (multiple sequence alignment)?
- A protein is a sequence over 20 amino acids (e.g. `MKTAYIAK…`). Over evolution, the same protein family has many "homologous sequences" across species: similar in function and structure, but mutated at certain positions.
- An MSA **aligns these homologous sequences position by position** and stacks them into a matrix:
  - **Row** = one homologous sequence; there are typically tens to tens of thousands of rows. Row 1 is the "query sequence" whose structure is to be predicted.
  - **Column** = the same residue position after alignment, i.e. the same "amino-acid slot". Positions that cannot be aligned are filled with gaps (`-`).
- Two kinds of signal can be read from this matrix:
  1. **Conservation**: if a column barely changes across species, that position is likely critical for function or structure.
  2. **Co-evolution / covariation**: if columns i and j always "change together" (when i mutates, j mutates as well), the two residues are often in contact in 3D. Early DCA and contact-prediction methods were built on statistics of this column-to-column covariation.
- Key point: **the information in an MSA is not in any single row but in comparisons across rows.** Looking at one sequence alone, you cannot tell which positions are conserved or which co-vary.

### Row / column (axial) attention in MSA Transformer
- The input is an S×L token matrix (S sequences, L positions each), with a d-dimensional vector per token.
- Full self-attention over all S·L tokens costs O((SL)²), which is too expensive. Hence an **axial factorisation**: each layer alternates between two kinds of attention:
  - **Row attention**: within each row (one sequence), the L positions attend to each other, i.e. ordinary sequence attention. MSA Transformer **sums the attention logits across all rows (tied row attention)**, yielding a single L×L attention map shared by all sequences; this map is itself approximately a "contact map".
  - **Column attention**: within each column (one position), the S sequences attend to each other. This is the "cross-sequence comparison": at position j, the query sequence can see which amino acid other species place at the same position.
  - The costs are roughly O(S·L²) and O(L·S²) respectively, far smaller than O((SL)²).
- AlphaFold2/3 do this slightly differently: in addition to the MSA representation, they maintain a separate **pair representation P∈R^{L×L×c}** (one vector per residue pair):
  - **Outer-product mean**: for each position pair (i,j), take the outer product of each row's vectors at i and j, **average over all rows**, and add the result to P_ij. This is a kind of "learnable inter-column covariance", a neural version of co-evolution statistics.
  - **Pair-weighted averaging** (AF3's MSA module): P is projected into L×L attention weights that mix information across positions **within each row**. The attention weights come from the pair representation, not from Q·K. AF3's MSA module no longer contains column attention.
  - Afterwards, the Pairformer processes only the single-sequence (single) representation and the pair representation, using the pair as an attention bias.

### How CellMSA maps this to cells × genes
| Protein MSA | CellMSA |
|---|---|
| One row = one homologous sequence | One row = one cell (row 1 is the target cell; the other 40 rows are retrieved context cells) |
| One column = one aligned residue position | One column = one gene (all rows use the same set of genes, i.e. they are "aligned") |
| Token = amino acid | Token = gene ID embedding + expression-bin embedding (11 bins, with 0 in its own bin); context rows additionally get a relation-type embedding |
| Conserved position | Marker genes stably expressed within a cell type |
| Co-evolving residue pair | Gene pairs that co-vary across cells (co-expression) |
| Pair representation L×L | Gene-pair representation G×G×d_p (d_p=8; G up to 2,048 genes plus 1 cls) |
| Homologous-sequence search (HHblits/jackhmmer) | Metadata-based retrieval: 16 of the same type and same batch, 16 of the same type but a different batch, 8 of related types |

**An important clarification:** CellMSA implements an **AF3-style MSA module** (outer-product mean + pair-weighted averaging + a Pairformer that only sees the target cell), **not** MSA-Transformer-style row/column axial attention. In the code (`cellmsa/model/msa_module.py`) each layer consists only of `OuterProductMean → MSAPairWeightedAveraging → Transition`; the triangle updates have been removed, leaving only a pairwise-transition linear layer. Therefore:
- There is **no attention at all** between cells (rows). Cross-cell information enters the pair representation **only** through "outer products averaged over rows". This is a set summary that is invariant to row order (permutation-invariant).
- The target cell ultimately "sees" the context only through a gene×gene bias matrix.
- This is a **different paradigm** from the "concatenated cells + bidirectional cross-cell attention" design described above. CellMSA can serve as a low-cost alternative or complement: compress the context into a gene-pair prior rather than letting cells attend directly to one another.

---

## A. CellMSA: Context Modeling for Single-Cell Representation Learning

- arXiv:2609.38908v1 [q-bio.GN], submitted 2026-09-30; authors Suyuan Zhao, Minghao Liu (co-first authors), Yizhen Luo, Zaiqing Nie (corresponding author); affiliations Tsinghua AIR and PharMolix. The Comments field on the arXiv page reads "Accepted by NeurIPS 2026, code released"; the article licence is CC BY 4.0.

### A.1 Problem setting and data
- **Task scope**: general-purpose single-cell representation learning. Downstream evaluation covers batch integration, cell-type and cell-state classification, and perturbation prediction (the latter via STATE-ST).
- **Pre-training corpus**: from CELLxGENE, about **109 million cell observations**, of which **65.6 million are primary observations** (`is_primary_data`); drawn from **2,090 datasets and 818 cell types**. The paper states explicitly that 109 million counts observations rather than unique cells; non-primary observations come mainly from re-aggregation within a study or cross-study atlas integration. Downstream datasets were removed from the pre-training corpus. Gene vocabulary: 61,982.
- **Downstream data**:
  - Batch integration: Tabula Sapiens, 483,152 cells, 24 tissues, 15 donors, 475 cell types; evaluated per tissue and macro-averaged.
  - Classification: Tabula Sapiens-Blood (cell-type annotation); PT cells from the Kidney Atlas, with three states aPT/dPT/dPT-DTL (77 donors in total: 26 healthy, 14 AKI, 37 CKD). Split 7:1:2 by donor, repeated 5 times.
  - **Perturbation prediction**: the Replogle dataset, **4 human cell lines** (K562, Jurkat, RPE1, HepG2). Only "the top 100 perturbations present in all 4 cell lines with the most cells" are kept; together with controls this gives about **132k cells**.
- **Perturbation-task split (important)**: in HepG2, **45% of perturbations are used for testing and 5% for validation**; **the remaining 50% of HepG2 perturbations + all perturbations from K562/Jurkat/RPE1 are used for training**. Because these 100 perturbations exist in all 4 cell lines, every test perturbation has been seen in the other 3 cell lines. This is therefore a "target cell line partially visible + test perturbations visible in other cell lines" setting, **neither zero-shot in cell line nor zero-shot in perturbation**.
- Cells per set: the paper does not state the STATE-ST cell-set length used in the perturbation task (unverified). The context pairs each target cell with 40 neighbours.

### A.2 Input construction (how context is selected and arranged)
- The context of each target cell c_i is 𝒩_i = 𝒩^same batch ∪ 𝒩^cross batch ∪ 𝒩^cross type:
  - **same batch**: same cell type, same batch; used for local denoising;
  - **cross batch**: same cell type, different batch; used to extract patterns that remain stable across batches;
  - **cross type**: biologically similar but different cell types, used as a "contrastive background" to highlight what is distinctive about the target cell.
- **Counts** (pre-training): 16 same-batch, 16 cross-batch, 8 cross-type, **40 neighbours** in total, i.e. S = 41 rows (row 1 is the target cell). In the code, `demo/prepare_context.py` uses the same widths of 16/16/8 for `stsb/stdb/dt`; at inference `donor_id` serves as the batch.
- **How related types are defined**: for each cell type, compute a cell-count-weighted mean across all batches (after normalising to 1e4 and log1p), compute cosine similarity across 61,982 genes, run Ward hierarchical clustering, and cut into **28 clusters** (k≈√818); within its cluster, each type takes up to **50** related types ranked by similarity. Cluster sizes range from 2 to 144; mean within-cluster cosine similarity is 0.845, between-cluster 0.581.
- **Relation-type embedding**: each neighbour carries a label r ∈ {same batch, cross batch, cross type} with a corresponding learnable E_rel, added to every token in that row; the target row gets none.
- **Label-based vs label-free retrieval**: the main batch-integration experiments retrieve using cell-type labels (as does Stack). For the classification task, the validation and test sets **do not use labels** and instead retrieve via HVG→PCA→KNN to avoid label leakage.
- **Context definition specific to the perturbation task**: context is drawn only from the same cell line. **Cells with the same perturbation** are treated as same-batch context, and **cells from the same cell line with a different perturbation** as cross-type context; **post-perturbation expression of test and validation perturbations may not be used as context**. The paper does not say how the cross-batch slots are filled in the perturbation task (unverified).
- **Tokenisation** (following scGPT): gene ID plus expression bin; 0 has its own bin and non-zero values are binned by within-cell quantiles, **11 bins** in total; maximum sequence length 2,048; a `<cls>` token is prepended.
- **How the gene axis is aligned (from the code, `cellmsa/dataset.py: select_gene_indices`)**: for each **target cell**, select its **most highly expressed non-zero genes**, up to 2,048 (`gene_truncation_nonzero_ratio=1.0`, so genes that "differ most from the neighbours" are not used as filler; any shortfall is simply padded), sorted by gene index; **context rows use the same set of gene columns**. In other words, the column set is determined by the target cell, and context-cell expression outside these genes is discarded.
  - (our analysis) For CRISPRi, this means that a target gene knocked down to 0, as well as genes that are 0 in the target cell but expressed in the context, **will not appear in the token sequence at all**. For weak-effect prediction, this "columns determined by non-zero genes" design should not be adopted as-is.

### A.3 Architecture details
**Overall**: input embedding → CellMSA-Module (low-dimensional, processes all S rows, outputs a gene-pair representation) → GenePairformer (high-dimensional, encodes only the target cell, pair used as attention bias) → `<cls>` taken as the cell embedding z_i∈R^512.

- **Input**: e_sg = E_gene(g) + E_value(x_sg) (+ E_rel(r_s), only for s>1).
- **Initial pair**: p⁰_uv = f_pair(E_gene(u), E_gene(v)). In code this is `Linear_left(e_u) + Linear_right(e_v)` followed by one Linear layer; pair dimension **d_p = 8**; the pair tensor has shape (b, G+1, G+1, 8), including cls.
- **CellMSA-Module** (L_MSA = **4** layers in total; e is first projected to a low dimension **d′ = 128**, the "MSA hidden size"). Each layer:
  1. **Outer-product Mean**: p^{l+1}_uv = p^l_uv + W_p · (1/S) Σ_s a_su ⊗ b_sv, where a and b are linear projections of m with projection dimension **8** (the "outer product hidden size"). **This is the only place where cross-cell interaction happens**, and it is an average over rows.
  2. **Pair-weighted Averaging**: A^h_uv = softmax_v(W_b^h p_uv), m_su ← m_su + W_o Concat_h(γ_su ⊙ Σ_v A^h_uv v_sv); **8 heads, 32 dimensions each**; γ is a sigmoid gate; row dropout 0.15. The attention weights come entirely from the pair representation, **all rows share the same gene×gene attention map**, and there is no Q·K.
  3. Transition (SwiGLU).
  - In the code there is a learnable layerscale at the end, and one pairwise-transition linear layer is kept; triangle multiplication and triangle attention are both **removed**.
- **GenePairformer** (**6** layers; hidden 512; pair-biased attention with **8 heads of 32 dimensions each**; row dropout 0.25):
  - R^{L+1,H}_uv = R^{L,H}_uv + Q_u K_v^T/√d_H (Uni-Mol style: logits are accumulated back into the pair);
  - Attention = softmax_v(R^{L+1,H}_uv) V_v.
  - Only the target cell's G gene tokens are processed.
- **Scale**: **47.13M** parameters in total; hidden 512; pair embedding 8; cell embedding 512.
- **Complexity** (paper Appendix A.2):
  - Outer-product mean: O(S·G²·d_p²); pair-weighted averaging: O(G²·d_p + S·G²·d′); GenePairformer per layer: O(G·d² + G²·d).
  - Total additional overhead O(L_MSA(S·G²·d_p² + S·G²·d′)), **linear in the number of context rows S** and quadratic in the number of genes.
  - (our analysis) Taking G=2,049 and S=41, the pair tensor for a single sample has about 4.2M×8 elements; pair-weighted averaging has to process S·G² ≈ 170 million elements times d′, so G² is the real bottleneck. Using a full transcriptome of about 18k genes would make this design infeasible; the gene count would first have to be reduced (HVGs or gene modules).

### A.4 Training objective and procedure; inference
- **L = L_MGM + λ₁·L_Rec + λ₂·L_CCE**, with λ₁ = **0.01** and λ₂ = **10**.
  - **MGM** (masked gene modelling): mask probability 0.15 (of which 80% replaced by mask and 10% randomly replaced); cross-entropy on the original bin. **The target row and the context rows are masked at the same gene positions**, so the model cannot simply copy the answer from neighbours.
  - **Rec**: reconstruct the bins of a random subset of genes from z_i plus E_gene(g) (reconstruction probability 0.1).
  - **CCE**: InfoNCE; positives are cells of the same cell type, negatives are other cells from the same batch; temperature τ = 0.05. The paper acknowledges that CCE uses cell-type metadata as weak supervision.
  - Appendix B.3 gives only a qualitative UMAP comparison: when CCE dominates, disease-state information is erased; when Rec dominates, the embedding is too noisy.
- **Optimisation**: AdamW, linear warmup for 1,000 steps then constant 1e-5, weight decay 0.05, mixed precision. Table A.32 in the paper lists a per-GPU batch of 4, 4 GPUs, gradient accumulation 2, i.e. an **effective batch of 32**; the code's `config.py` has batch_size 16 and gradient_accumulation 16. These are inconsistent (unverified which is the final setting).
- **Compute**: **4×A800 for 1 epoch, about 20 days**.
- **Inference**: given a query h5ad, a reference h5ad, and context indices for each query cell (`stsb/stdb/dt` or `neighbor_indices`) → output the cls embedding. Downstream tasks train an MLP on the embedding or attach STATE-ST.
- **Perturbation-prediction pipeline**: CellMSA is frozen as an encoder to produce cell representations; **STATE-ST** is then trained with the official STATE implementation and evaluated with **Cell-Eval**. Because STATE-ST must be retrained each time, which is costly, **only a single fixed split and a single run were performed** (no variance estimate). The paper does not state whether Cell-Eval is computed in gene space or embedding space, nor how embeddings are decoded back to expression (unverified).

### A.5 Results
**Perturbation prediction (Table 3, Replogle subset, HepG2 test)**
| Model | Pearson Δ | PRAUC | Spearman-FC | DE Overlap |
|---|---|---|---|---|
| Expression + STATE-ST | 0.398 | 0.276 | 0.404 | 0.139 |
| STATE-SE + STATE-ST | 0.353 | 0.287 | 0.365 | 0.180 |
| Stack + STATE-ST | 0.358 | 0.293 | 0.372 | 0.148 |
| **CellMSA + STATE-ST** | **0.433** | **0.334** | **0.431** | **0.215** |

- Metric definitions (Appendix C.5, Cell-Eval): Δ_p = perturbed-group mean − control-group mean, and Pearson Δ is the correlation across genes; PRAUC uses true DE significance as labels and ranks by predicted DE significance; DE Overlap takes the top DE genes by |FC| and computes the overlap; Spearman-FC is computed only on truly significant DE genes.
- The "8.8% improvement over the strongest baseline" = 0.433 / 0.398 (the strongest baseline is **raw expression**, not STATE-SE or Stack). Notably, **both foundation-model embeddings (STATE-SE, Stack) underperform raw expression on Pearson Δ**.

**Ablations (Table 4)**
| Variant | PT classification Acc | PT Macro F1 | Perturbation Pearson Δ | Perturbation Spearman-FC |
|---|---|---|---|---|
| w/o CellMSA-Module (reduces to a plain Transformer) | 0.854 | 0.811 | 0.355 | 0.366 |
| w/o Context (same architecture, context rows removed) | 0.878 | 0.817 | 0.389 | 0.367 |
| Only Same-batch Context | 0.929 | 0.893 | 0.400 | 0.397 |
| **Full** | **0.958** | **0.931** | **0.433** | **0.431** |

- (our analysis) On the perturbation task, "w/o Context" (0.389) is slightly below the raw-expression baseline (0.398); only the variants with context exceed raw expression. The increment from cross-batch and related-type context (0.400→0.433) is larger than the increment from having context at all (0.389→0.400). This is, however, a single run without error bars.
- **Context size (Figure 5, PT classification only)**: performance first rises and then saturates, with "no substantial improvement beyond 40 context cells". The exact values in the figure are not given in text form (unverified).

**Other tasks (brief)**
- Batch integration (Table 1, macro average): CellMSA Total **0.736**, Stack 0.663, scVI 0.600, STATE-SE 0.592, scGPT 0.586, Geneformer 0.584, PC-HVG 0.578; CellMSA's Bio score is 0.850 and Batch score 0.566. CellMSA and Stack use label-based retrieval and thus have access to different information than the other methods.
- Label-free retrieval (Bladder, Table A.27): CellMSA Total **0.688**, Stack 0.600, Geneformer 0.612, scGPT 0.598; switching to label-based retrieval raises CellMSA to 0.748.
- Classification (Table 2): cell type Acc 0.962±0.002, Macro F1 0.912±0.005 (CellPLM: 0.931 / 0.906); PT state Acc 0.958±0.004, Macro F1 0.931±0.006.
- Robustness: after randomly dropping 10%, 30%, 50% of non-zero genes, the mean embedding Pearson is 0.996, 0.988, 0.972 respectively.
- Interpretability: heads 6/7 favour disease-related gene pairs, heads 4/5 favour homeostasis-related ones; 26 of KDR's top-50 paired genes (52%) are recorded in STRING v12. The authors note this is descriptive only, with no significance test.

**Baselines**: PC-HVG, scVI, scGPT (whole-human), Geneformer V2-104M, STATE SE-600M, Stack-Large, CellPLM, all used as zero-shot embedding extractors without fine-tuning. **The perturbation task is not compared against dedicated perturbation models such as GEARS, CPA, or linear baselines.**

### A.6 Code and weights
- GitHub: https://github.com/PharMolix/CellMSA (created 2026-09-29). **The LICENSE is MIT**, noting that the third-party components alphafold3-pytorch and hyper-connections are copyright Phil Wang (the GitHub API reports "Other/NOASSERTION", but the LICENSE file body is the MIT text).
- Contents: `best_model.pt` (stored with Git LFS, containing the weights and gene vocabulary; actual file size unverified), `embedding.py/.sh` (inference), `demo/` (5,000 blood cells + `prepare_context.py`), `mlp_classifier.py`, `cellmsa/model/{cellmsa.py, msa_module.py, attention.py}`, `dataset.py` (including the PrebuiltDataset used for pre-training).
- **No pre-training entry script and no perturbation-prediction (STATE-ST) pipeline were found in the repository tree** (as of 2026-10-04). The code's `build_neighbor_relation_ids` splits relation IDs as 24/24, inconsistent with the paper's 16/16/8; the inference path, however, uses 16/16/8 (which was actually used during training is unverified).
- Environment: Python 3.12, PyTorch 2.5.x; inference requires an NVIDIA GPU.

### A.7 Implications for the in-context design (CellMSA)
**Worth adopting**
1. **Row-role embeddings**: add a "relation-type" embedding to each context cell. For the design under consideration, possible roles are {new-cell-line control; training-cell-line control; training-cell-line cell perturbed by target t′; …}, combined with "similarity bucket between t′ and the query target t" and "the cell line t′ comes from". This is far more reliable than relying on positional encodings alone, and lets the model learn to distinguish "this row is background, that row is a reference response".
2. **Context should be diverse, not necessarily large**: going from same-batch-only context to adding cross-batch and related-type context gives the largest increment; saturation occurs at about 40 cells. Rather than packing in thousands of cells, it is better to curate: a number of new-cell-line controls + a small number of perturbed cells for pathway-related targets across multiple cell lines, with diverse sources.
3. **Joint masking at the same positions to prevent copying**: CellMSA masks the same genes in the target and the context. The corresponding risk here: if, during training, the context includes perturbed cells from "the same cell line + the same target", the model will learn to copy the answer, yet such context is unavailable at test time. **Training construction must strictly mimic test-time visibility**: no perturbed cell from the query cell line may appear in the context, and sampling should be leave-cell-line-out.
4. **Summarise context cheaply with an outer-product mean**: if full cell×cell attention (multiplied by gene tokens) is too expensive, one can take the outer-product mean of "Δ expression for t′ in other cell lines" over the context, obtaining a gene-gene "response-coupling prior", and inject it as a bias into the gene encoder of the new cell line's control cells. The cost is linear in the number of context rows and independent of order.
5. **Two-stage asymmetric structure**: first process many context rows cheaply in a low dimension, then process the target in detail in a high dimension. Here this could mean: "lightweight encoding of many context cells → heavy modelling of the new cell line's controls and their output distribution".

**Not directly transferable**
1. CellMSA's perturbation experiment: (1) 50% of the target cell line's perturbations are visible; (2) test perturbations have been seen in all 3 other cell lines; (3) only the 100 perturbations with the most cells are used (biased toward strong effects); (4) single split, no error bars. **There is no evidence that it works in a zero-shot cell-line, weak-effect setting.**
2. Its perturbation context is **drawn only from the same cell line**. Our setting is the opposite: the new cell line has no perturbed cells at all, so perturbed cells can only be borrowed from other cell lines. CellMSA does not study "borrowing perturbation responses across cell lines".
3. **The gene axis is determined by the target cell's non-zero genes**: for CRISPRi, a target gene that is knocked down may vanish from the tokens entirely; and weak-effect signals often lie in low-to-moderately expressed genes, which the 2,048 truncation discards. A **fixed, aligned gene panel** is needed, and it must include the target genes themselves.
4. The output is an embedding (cls), not an expression distribution. Generating a population of perturbed cells requires a separate decoder (NB/ZINB, flow, etc.).
5. Pre-training used 109 million observations and 20 days on 4×A800; the MIT-licensed weights can be used as initialisation or features, but its tokenisation (non-zero truncation, 11 bins) is lossy for count-level modelling of subtle changes.

**Open questions**
- The outer-product mean summarises **covariance between cells**. For weak-effect perturbations, Δ is much smaller than cell-to-cell heterogeneity, so the summary is likely to be dominated by cell-state variation. Should one pseudo-bulk first, or take differences (perturbed cell − matched control) before summarising?
- Context rows can interact with the pair representation only via a "mean over rows", so the model cannot express selectivity of the form "only one particular row is relevant". This is exactly the problem PT-RAG tries to address with Gumbel selection.

---

## B. Retrieval-Augmented Generation for Predicting Cellular Responses to Gene Perturbation (PT-RAG)

- arXiv:2603.07233v1 [cs.LG], submitted 2026-03-07; authors Andrea Giuseppe Di Francesco (Sapienza/ISTI-CNR), Andrea Rubbi (Cambridge/Sanger; the two are co-first authors), Pietro Liò (Cambridge). Comments: "Accepted at ICLR 2026 Workshop: Generative AI in Genomics. 25 pages, 9 figures"; article licence CC BY 4.0.

### B.1 Problem setting and data
- **Task**: given a population of control cells and a perturbation ID, predict the distribution of perturbed cells (single-gene perturbations).
- **Data**: Replogle-Nadig (Replogle 2022 + Nadig 2024), 4 cell lines: **K562, Jurkat, RPE1, HepG2**. **2,009 unique perturbations** in total (intersection of Replogle-Nadig with genes that have a GenePT embedding), **2,000 HVGs**.
- **Split ("few-shot cross-cell-type")**: for each target cell line, train on the other 3 cell lines and add **30% of the target cell line's perturbations** to training as few-shot examples; the remaining **70%** are used for validation and testing ("similar to Adduri et al. 2025"). Final results are averaged over the 4 cell lines.
- **Test size**: **1,635 test perturbations** in total (HepG2 375, RPE1 416, Jurkat 443, K562 401).
- **Cells per set**: STATE's cell set is **64 cells** (`cell_set_len=64`), batch size 64.
- **Conclusion**: **not a zero-shot cell-line setting** (30% of the target cell line's perturbations are visible). Whether the test perturbations appear in other cell lines is not stated case by case in the paper; given "training on all perturbations of the other 3 cell lines", they very likely do (unverified).

### B.2 Input construction and retrieval
- **Perturbation representation**: instead of one-hot, **GenePT** is used: GPT-3.5 (text-embedding) encodings of NCBI gene descriptions, **1,536 dimensions**, normalised before computing cosine similarity.
- **Retrieval corpus**: the **GenePT embeddings of all 2,009 perturbations** in the training set. The query perturbation itself is excluded at retrieval (in code, its self-similarity is set to −inf).
- **(Key fact) The retrieved "content" is only the GenePT vectors of the candidate perturbations**, i.e. h^cxt_k = PertEncoder(h^gpt_{p(k)}); **it does not include the measured expression responses of those perturbations in any cell line**. The introduction speaks of "conditioning the generator on observed responses to related perturbations", but in the method section (Eq. 4) and the code (`refine_combined_input_diff_retrieve_than_predict`) the retrieved items are always embeddings. So the "RAG" here is in practice **a learnable gated mixture of the embeddings of GenePT-neighbour genes**.
- **Two stages**:
  1. **Semantic retrieval**: top-**K=32** by GenePT cosine similarity (no similarity threshold; always 32), non-differentiable.
  2. **Differentiable selection**: for each candidate, build a triple c_k = [h^ctrl ; h_pert ; h^cxt_k] (3×128=384 dimensions) → LayerNorm → MLP_score (2 layers, hidden 128) → 2 logits (exclude/include). In the code the logits are further normalised by the std of the whole tensor → **Straight-Through Gumbel-Softmax** (τ = **0.5**, hard=True), yielding w_k∈{0,1}.
- **Cell-type awareness**: what the authors call "cell-type awareness" refers solely to the fact that the scorer's input includes h^ctrl (the embedding of **each control cell**). In the code, selection happens over the (B, S, K) dimensions, so **each of the 64 cells independently picks its own context subset**. The candidate pool (GenePT top-32) is independent of cell type.

### B.3 Architecture details
- **Backbone**: STATE (State Transition). The paper says "Llama backbone, sequence (cell population) 64 cells, batch 64".
  - The repository's default config `configs/model/state.yaml`: 8 layers, 12 heads, `bidirectional_attention: false`, hidden 696; but the experiment script `run_celltype_experiments.sh` overrides this to **`hidden_dim=128`**, `cell_set_len=64`, `batch_encoder=true` (`gem_group` as batch), `predict_residual: True`, `loss: energy`, `embed_key=X_hvg`, `output_space=gene`.
  - (unverified) Whether the experiments override `bidirectional_attention`. If not, this Llama backbone applies **causal (unidirectional) attention** across the 64 cells rather than bidirectional attention. This is directly relevant to a bidirectional cross-cell-attention design and is worth double-checking.
- **Cell Encoder**: according to the paper, the pre-trained SE-600M (HuggingFace `arcinstitute/SE-600M`) maps 2,000 genes to **128 dimensions** and is **frozen throughout training**. The script uses `init_from=se600m_epoch16.ckpt` plus `freeze_pert_backbone=true`, together with `embed_key=X_hvg`; exactly how SE-600M is wired to the HVG input could not be fully determined from the code (unverified).
- **Perturbation Encoder**: single-layer MLP, 1,536 → 128.
- **Generation**:
  - Baseline: z = h^ctrl + h^pert → Transformer → x̂^pert.
  - Vanilla RAG: z = CrossAttention(q = h^ctrl + h_pert, k = v = h_cxt). In code: `combined = query_enc(combined)`, then `combined_input = cross_attention(query=combined, key=…, value=…)[0]`, **with no residual, directly replacing** z. That is, **the output z is merely a convex combination of the context values**; the control-cell state and the query perturbation influence it only indirectly through the attention weights.
  - PT-RAG: h′_k = MLP_proj(c_k) (single layer, hidden 128), z = Σ_k w_k h′_k, which **likewise directly replaces** h^ctrl + h^pert. Because the triple contains h^ctrl and h_pert, that information can be preserved through MLP_proj; but if a cell selects no candidate at all, z = 0 (our analysis, inferred from the code).
- **Parameters and FLOPs (Table 7)**: STATE 20,936,480 / 1.67B FLOPs per batch / 34.9M FLOPs per cell; STATE+GenePT 20,874,016 / 1.67B / 34.8M; Vanilla RAG 21,120,288 / 1.67B / 34.8M; **PT-RAG 20,973,602 / 2.86B / 59.9M** (about 1.7×).
- **Complexity**: each cell scores K candidates, O(B·S·K), plus the STATE backbone's O(S²) attention over S cells.

### B.4 Training objective and procedure; inference
- **L = L_dist + λ_sparse·L_sparse**: L_dist is the **energy distance** (inherited from STATE); L_sparse = (1/K) Σ_k w_k (L1 sparsity), **λ_sparse = 0.1**.
- **Optimisation**: Adam, lr **1e-3**, weight decay **0.0005**, up to **50,000 steps**, validation every 2,000 steps; typically converges in 30,000–40,000 steps.
- **Compute**: A100 40GB, **about 8–10 hours per target cell line**.
- **Inference**: Gumbel uses hard argmax in the forward pass. The code also calls `F.gumbel_softmax(..., hard=True)` at inference, so **Gumbel noise is present at inference too and selection is stochastic** (our reading of the code; the paper says "maintaining hard selections at inference time"). The checkpoint and retrieval index are saved, and inference is run with `state tx predict`.

### B.5 Results
**Main table (Tables 1/2, averaged over 4 cell lines, n = 1,635; std in parentheses)**
| Metric | STATE | STATE+GenePT | Vanilla RAG | PT-RAG |
|---|---|---|---|---|
| Pearson DEG ↑ | 0.624 (±0.048)† | 0.631 (±0.051) | 0.396 (±0.063)† | **0.633** (±0.048) |
| Spearman DEG ↑ | 0.403† | 0.411 | 0.307† | **0.412** |
| MSE ↓ | 0.211 | **0.210** | 0.316† | **0.210** |
| RMSE ↓ | 0.458 | 0.458 | 0.562† | **0.457** |
| MAE ↓ | 0.298† | 0.296 | 0.429† | **0.295** |
| MSE_PCA50 ↓ | 8.43 | 8.42 | 12.64† | **8.39** |
| W1 ↓ | 35.70† | 35.53††† | 48.48† | **35.41** |
| W2 ↓ | 646.1† | 638.7†† | 1189.5† | **633.7** |
| Energy ↓ | 9.41††† | 9.40 | 14.18† | **9.33** |
| MMD ↓ | 0.0142 | 0.0142 | 0.0142 | 0.0142 |

(† pFDR<0.01, †† <0.05, ††† <0.1; test: Mann-Whitney U + BH-FDR, 3 pairwise comparisons × 10 metrics = 30 tests)

- **PT-RAG vs STATE+GenePT (Table 8)**: only **W2 is significant (pFDR = 0.041)**, with W1 marginal (0.082); Pearson DEG (0.590), Spearman (0.785), MSE/MAE, etc. are **not significant**. The authors themselves acknowledge that the gain comes mainly from the GenePT representation.
- **PT-RAG vs STATE**: Pearson DEG pFDR = 2.44e-8, Spearman 4.89e-10, MAE 9.33e-5, W1 3.41e-6, W2 5.47e-8, Energy 0.082.
- **Per cell line (Tables 3–6)**: PT-RAG is best on HepG2 (Pearson DEG 0.604 vs STATE 0.596) and RPE1 (0.640 vs 0.621). On Jurkat, STATE+GenePT is best (0.652 vs PT-RAG 0.649). **On K562, STATE is best** (0.639 vs PT-RAG 0.633; W2 576.9 [GenePT] vs 577.1).
- **Metric definitions**: DEGs are called with a Welch t-test (p<0.05) between control and perturbed groups, and Pearson/Spearman between "predicted change and observed change" are computed on the DEGs. Distributional metrics are computed in a PCA-50 space fitted on the training set. **Pearson DEG here is not the same metric as Cell-Eval's Pearson Δ** and cannot be compared directly with CellMSA's 0.433.

**Ablations and sensitivity (all on HepG2 with K=32)**
- **λ_sparse**: with λ = 0, on average **31.949/32** candidates are selected, Pearson DEG is only **0.134**, Spearman **−0.025**, RMSE 0.567, MSE 0.322, MAE 0.362, W2 651.744, Energy 11.974. With λ ∈ {0.01, 0.1, 1.0}, on average **12.536 / 6.612 / 4.915** are selected, and Pearson DEG is stable at 0.594–0.604, Spearman at 0.386–0.401.
- **K = 16 vs 32**: Pearson DEG of 0.578–0.604 and 0.594–0.604 respectively; little difference.
- **Vanilla RAG vs K** (K ∈ {2, 5, 10, 32}): Pearson DEG rises from **0.293** at K=2 to **0.351** at K=32 (PT-RAG: 0.604); Spearman goes 0.220 → 0.170 → 0.189 → 0.289; RMSE 0.608 → 0.631 → … → 0.577; MSE 0.370 → 0.399 → … → 0.334.
- **Ablations not performed**: random retrieval vs GenePT retrieval, different candidate sources, retrieving measured expression responses, zero-shot cell line. **A "retrieval vs random" comparison, which is central to the question of whether retrieval helps, is absent from the paper.**

### B.6 Code and weights
- GitHub: https://github.com/difra100/PT-RAG_ICLR (created 2026-02-09, last push 2026-08-03). It is a fork of Arc's STATE code; the PT-RAG logic is in `src/state/tx/models/state_transition.py`. The repository provides data-download scripts and per-cell-line split TOMLs (`datasets/repogle_nadig*.toml`, which give the few-shot splits as explicit perturbation lists).
- **The repository has no LICENSE file** (the GitHub API returns license: None), which by default means all rights reserved (should be confirmed before any reuse). The underlying STATE code and SE-600M weights are under Arc Institute's licences; their specific terms were not checked this time (unverified).
- No released trained PT-RAG checkpoint was found (unverified).

### B.7 Critical assessment of the paper's claims (our analysis)
1. **The Jaccard evidence for "cell-type-aware retrieval" is indistinguishable from random selection.** The paper reports that, for 33 shared genes, the Jaccard index between the top-10 selected perturbations of different cell lines is 0.185–0.196, with a mean of 0.191. However, the candidate pool is the **same** GenePT top-32 for every cell line (independent of cell type). For two independent random draws of 10 elements from 32, the expected Jaccard is **≈0.191** (our simulation, 200,000 draws). The reported cross-cell-type Jaccard therefore matches the expectation under random selection from a shared 32-candidate pool, and on its own does not demonstrate that cell-type-specific selection has been learned. The example in Appendix E.1 (WARS selecting other tRNA synthetases) likewise reflects the fact that GenePT's top-32 already consists of the same gene family, and does not by itself speak to the second-stage selection.
2. **The degradation of Vanilla RAG may stem from the implementation rather than from retrieval being harmful.** Its cross-attention output **replaces** h^ctrl + h_pert without a residual connection, so the Transformer's input is reduced to "a convex combination of context GenePT vectors". PT-RAG's triple contains h^ctrl and h_pert, so that information is preserved, which may be the real reason for the gap between the two. Similarly, when all 32 candidates are selected at λ = 0, the **sum** scales the representation by roughly 32×, which could cause numerical or scale problems (the paper attributes this to "noise"). The paper's conclusion that "non-differentiable retrieval is harmful" has not been tested against the control of "Vanilla RAG with a residual connection".
3. **MMD is exactly identical across all models and all cell lines** (e.g. 0.0142±0.0079), which may indicate an issue in the evaluation; this metric is best disregarded.
4. **Some numbers in the main text differ from the tables.** §4.4 states "Pearson: 0.293 vs 0.624; Spearman: 0.220 vs 0.403", whereas Table 1 lists Vanilla RAG as 0.396/0.307; 0.293/0.220 correspond to the K=2 ablation values on HepG2.
5. **The effect size is small.** PT-RAG's Pearson DEG exceeds STATE's by only 0.009 (std ≈0.048); the significance comes from an unpaired rank test with n = 1,635. Per cell line, STATE is actually better on K562.

### B.8 Implications for the in-context design (PT-RAG)
**Worth adopting**
1. **Represent targets with GenePT or other gene priors**, and use this for candidate pre-screening (top-K). This step captures nearly all of the gain (STATE → STATE+GenePT: Pearson DEG 0.624 → 0.631). In our setting, a combination of GenePT, STRING, GO, or essentiality could be used to choose the candidate pool of "pathway-related targets".
2. **Straight-through Gumbel + sparsity regularisation to learn "which context to include"**: if there are many reference targets t′ in the context, one can score quadruples such as [new-cell-line control state; query target t; candidate target t′; summary of t′'s response in its source cell line] and then make hard selections. A sparsity or budget constraint is required (λ=0 collapses).
3. **Context injection must have a residual connection, with the gate initialised to 0** (the lesson from B.7-2): when doing bidirectional attention over concatenated cells, the model should also be able to fall back to the no-context baseline when "none of the context is useful". For example, add a learnable gate initialised to 0 on the attention output from context tokens, or include a drop-all-context training branch.
4. **Normalise for context size**: use mean or softmax weighting instead of a sum, so that scale does not drift with the number of context items.

**Not directly transferable**
1. PT-RAG **never retrieves measured perturbation-response data**. Placing perturbed cells for related targets from other cell lines into the context is intrinsically more informative, and PT-RAG's results neither support nor refute that approach.
2. Its setting is few-shot (30% of the target cell line's perturbations visible), and most perturbations are strong-effect knockdowns of essential genes, whereas our setting is zero-shot in cell line with many weak effects. **The only entry point for "cell-type awareness" in PT-RAG is h^ctrl**; in the zero-shot setting, the new cell line's controls are likewise the only available cell-line information, but the paper does not show that this entry point actually has an effect (see B.7-1).
3. The backbone has only 128 dimensions, 2,000 HVGs, and sets of 64 cells. Modelling count-level distributions (rather than PCA-space or HVG summaries) would need a separate count decoder.
4. The evaluation metrics (Pearson on Welch DEGs, W1/W2 in PCA-50 space) are specific to this paper and are less sensitive to weak, low-expression effects, so the numbers do not transfer directly to other evaluation protocols.

**Open questions**
- For weak-effect targets, is "the response of related targets in other cell lines" more informative than "the prior embedding of the query target itself"? This calls for a three-way **retrieval vs random vs none** comparison under a leave-cell-line-out setting (exactly the comparison missing from PT-RAG).
- Should selection happen at the **set level** (all predicted cells for a target sharing one set of reference targets) rather than per cell? PT-RAG's stochastic per-cell selection injects additional variance into the generated distribution, which is not necessarily beneficial for energy/W-type metrics.

---

## C. Synthesis: concrete recommendations for an in-context predictor using concatenated cells + bidirectional cross-cell attention

1. **Context recipe (following CellMSA's 16/16/8 idea)**: new-cell-line controls (role A, analogous to same-batch; for denoising and characterising the cell-line baseline) + controls from other cell lines (role B, analogous to cross-batch; lets the model see "cell-line differences") + perturbed cells for related targets t′ from other cell lines, together with their matched controls (role C, analogous to the cross-type "contrastive background"; provides response templates). Each row gets a role embedding, a cell-line embedding, and an identity or similarity embedding for t′. Start with about 40–64 cells in total and run scaling ablations (CellMSA saturates at about 40; STATE/PT-RAG use 64).
2. **Training sampling must strictly mimic testing**: no perturbed cell from the query cell line may enter the context (leave-cell-line-out); perturbed cells for the query target t from other cell lines **may** enter the context (they are also available at test time). In addition, include a sub-task that "provides only new-cell-line controls", to prevent the model from degenerating into copying the context.
3. **Cost of cross-cell interaction**: full bidirectional attention is quadratic along both the cell and gene axes. Options: (a) encode cells into vectors first and then do cell-level bidirectional attention (STATE-style); (b) following CellMSA, summarise role-C context (perturbed minus matched control) into a gene-pair bias via an outer-product mean and inject it into the gene encoder of role-A cells; (c) combine both.
4. **Special handling for weak effects**: both CellMSA (non-zero gene truncation, 11-bin discretisation) and PT-RAG (2,000 HVGs, PCA space) discard subtle changes in low-to-moderately expressed genes. A **fixed gene panel** that always includes the perturbed target genes is needed; train in counts space (NB-type likelihood), or design a dedicated loss on Δ.
5. **Essential controls**: none / random-target context / GenePT retrieval / data-driven retrieval (based on response similarity in other cell lines), each under the leave-cell-line-out setting. Also report **variance across multiple splits** (missing from both papers).
6. **Directly reusable assets**: in CellMSA's MIT-licensed code, `OuterProductMean` / `MSAPairWeightedAveraging` / `PairformerStack` are clean reference implementations (based on lucidrains' alphafold3-pytorch); PT-RAG's Gumbel-selection code is short and can be reimplemented from its description (the repository has no LICENSE, so it **should not be copied directly**).

---

## Sources
- CellMSA arXiv abs: https://arxiv.org/abs/2609.38908
- CellMSA HTML full text: https://arxiv.org/html/2609.38908 (v1, 2026-09-30)
- CellMSA code: https://github.com/PharMolix/CellMSA (README, LICENSE, `cellmsa/config.py`, `cellmsa/dataset.py`, `cellmsa/model/msa_module.py`, `cellmsa/model/cellmsa.py`, `demo/prepare_context.py`, `embedding.py`)
- PT-RAG arXiv abs: https://arxiv.org/abs/2603.07233
- PT-RAG HTML full text: https://arxiv.org/html/2603.07233 (v1, 2026-03-07); PDF: https://arxiv.org/pdf/2603.07233 (the HTML version was used as the reference; it was not compared page by page against the PDF)
- PT-RAG code: https://github.com/difra100/PT-RAG_ICLR (README, `src/state/tx/models/state_transition.py`, `src/state/configs/model/state.yaml`, `run_celltype_experiments.sh`, `datasets/repogle_nadig.toml`)
- The background section (MSA Transformer / AlphaFold) is based on the descriptions in CellMSA §2, §3.2, and Appendix A, plus general domain knowledge; the original MSA Transformer paper (Rao et al., ICML 2021) was not re-fetched for this note.
- Jaccard random baseline: Monte Carlo simulation (two random draws of 10 from 32, 200,000 repetitions, mean 0.1918).
