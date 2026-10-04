# Stack (Arc Institute, bioRxiv 2026): a technical deep dive, read with perturbation prediction in never-perturbed cellular contexts in mind

> **Subject:** Dong, Adduri, …, Roohani. *Stack: In-Context Learning of Single-Cell Biology*. bioRxiv 10.64898/2026.01.09.698608.
>
> **Motivating setting:** this review was written for a long-term research project on how genes and cells influence each other across diverse cellular contexts. Its core evaluation problem is **predicting perturbation responses in cellular contexts never seen perturbed**, e.g. CRISPRi knockdown responses in a cell line for which only unperturbed control cells are available. The review asks what Stack does, how it does it, and what can be reused for that problem.
>
> **Basis:** mainly the **v1 (2026-01-09)** full text (PMC12803207), the supplementary PDF (from the Europe PMC supplementary-files API), the code repository `ArcInstitute/stack` at commit `cacc2e4` (2026-04-27) and the Hugging Face model cards.
>
> **Note: bioRxiv already has a v2 (2026-06-08, with the author list grown to 17).** Its abstract says Perturb Sapiens has been extended to "892 drug, cytokine, and genetic perturbations" and adds a DiseasePert-3M dataset. The v2 full text was behind a Cloudflare/proof-of-work wall and could not be read, so **every v2 change to methods or results is (unverified)**. All numbers below come from v1.
>
> **Code references** use the form `path:line`, with paths relative to `src/stack/` in the `ArcInstitute/stack` repository at commit `cacc2e4`.
>
> **Markers:** **[code]** = confirmed by reading the code; **[paper]** = stated in the paper; **[inference]** = the author's own reasoning; **(unverified)** = could not be confirmed.

---

## 0. TL;DR (for decision makers)

1. Stack = "project each cell into 100 gene-module tokens (d=16)" + "tabular block: MHA across the tokens within a cell → MHA across cells (with the 100×16 = 1,600-dim flattened cell treated as one token) → token-wise FFN" × 9 layers + "a cell-wise MLP that decodes to a negative binomial (NB) distribution (scVI-style softmax × library size)". Large has 217M parameters, of which ~101M sit in the decoder and ~92M in the 1,600-dim projections of the 9 inter-cell attention layers. A hand count of the parameters from the code **matches Table 3 of the paper exactly, entry by entry** (G = 15,012).
2. Pre-training: 148.8M cells from scBaseCount. Each SRX file is cut into **contiguous 256-cell blocks**, and each block is one cell set. A rectangular gene mask is used (the same genes are masked in every cell at once, mask ratio ~U(0.1, 0.8)). The loss is NB NLL + 0.01 × sliced Wasserstein (pulling the set-centred embeddings towards N(0, I)).
3. Post-training ("Aligned"): the cell set of one sample is sorted by cell type. The first 25% is the **prompt condition** (causal mask: these cells can only see each other). Of the remaining 75%, part is "prompt context" (ground-truth cells, fraction ~U(0, 0.9375)), and the remaining positions are **replaced by cells of the same cell type from another sample (the query)**. The student model predicts the original cells (the target) at the replaced positions. Loss = energy distance in gene space (log-CP10k, top 1,000 Pearson-residual HVGs, stratified by cell type) + energy distance in embedding space (aligned to the target embeddings of an EMA teacher) + reconstruction NLL + SW + a classifier that tells prompt context from query (used for confidence at inference time). In essence this is **masked diffusion with cells as tokens**.
4. Inference: `[prompt cells 25% | prompt-context cells 20%→40% | query cells]` are concatenated into a 512-cell window (the window size is an inference). Over T = 5 steps the query is iteratively "unmasked", and the final output is **integer counts for every query cell (sampled from the NB)**, on a fixed 15,012-gene axis.
5. Key facts for predicting perturbations in never-perturbed contexts:
   - **The v1 post-training data contain no genetic perturbations at all** (45M observational cells from CELLxGENE + 10M Parse PBMC cytokine cells). The paper's only genetic-perturbation result is an **embedding linear probe** on X-Atlas/Orion (Xaira) (classification of 50 perturbations). Its balanced accuracy is ~0.04 (chance = 0.02), and only one setting, Stack-Large (All Context), beats the PC-HVG baseline. **There is no ICL generation or prediction result for genetic perturbations** (v1).
   - The code already has a `type="drug"` dataset config: group = condition, identity = cell line, query = **controls from the same cell line**. That is exactly a cross-cell-line transfer episode ("the same perturbation in other cell lines → controls of a new cell line"), and it **could be used directly to post-train on CRISPRi data**. However, it requires the prompt and the target to share **the same condition**. It does not support a prompt made of "pathway-related targets", and it has no perturbation-identity input.
   - Gene vocabulary: Stack's vocabulary is a fixed set of **15,012** genes, chosen as the union of the per-file top-1,000 HVGs. Before HVG selection, a gene-name regex filter removes many genes; one pattern, `.*P[0-9]*$`, wrongly drops ordinary protein-coding genes such as **TP53, MMP9, TIMP1, IGFBP3, AQP1, BMP4, GBP1**. Any gene outside the vocabulary is invisible to the model in both input and output; if a knocked-down target falls outside it, the knockdown of the target gene itself cannot be seen or predicted.
   - Licence: the code is CC BY-NC-SA 4.0. The weights and **outputs** are under Arc's non-commercial licence, and "Derivative Work" explicitly includes fine-tuning, distillation and **"use as a prior"**. Commercial or company-affiliated use requires confirmation against the licence **(unverified)**. Nothing in this document is legal advice.
6. Recommendation: **do not use Stack directly as the main prediction model.** The most useful things to borrow are (a) the tabular block design (~130 lines to re-implement inside a flow-matching model, swapping in SDPA/flash attention); (b) the episode construction of "cell-as-token masked diffusion + causal prompt block + query position embedding"; and (c) the `drug`-type cross-cell-line episode sampler. If the weights are used at all, limit it to offline comparison experiments (frozen encoder / synthetic-control Δ), and clarify the licence first.

---

## 1. Data and preprocessing

### 1.1 Pre-training corpus
- [paper] scBaseCount (Youngblut et al. 2025, human data uniformly re-aligned from SRA) has 189M cells after strict QC. The training set is **19,978 SRX samples, 148.8M cells** (Table 2; the main text says 149M), with the remaining ~20% held out for validation/test. A CELLxGENE version (905 datasets, 73.7M cells) and an scBaseCount-subset version (9,004 samples, 60.2M cells) were also trained.
- [paper] Cell filtering: scBaseCount keeps cells with 300–7,000 detected genes and ≥700 UMIs. CELLxGENE is not filtered.
- [code] The pre-training split is **by file** (`data/training/datasets.py:537-569`, `_split_files`, with test/val defaulting to 0.1 each via the CLI), not by cell.
- Whether scBaseCount contains the SRX runs of Perturb-seq datasets such as Replogle or X-Atlas is **(unverified)**. This matters for whether the cell lines used for evaluation have already been "seen".

### 1.2 Gene axis / vocabulary
- [paper] The unified gene list = for each scBaseCount file, the top 1,000 HVGs by **analytic Pearson residuals**, taken as a union and capped at 15,012.
- [code] `data/hvg.py:34-48` computes analytic Pearson residuals (θ = 100; residuals are not clipped, only NaN/inf are set to 0). `hvg.py:160-171` takes the 1,000 genes with the largest variance per file and stops once the union exceeds 15,000. By default `filter_gene_names_flag=True` (`hvg.py:58`) applies the regexes in `data/gene_processing.py` to remove: `^ENS…`, `^LOC…`, `^LINC…`, `^MT-`, `^RP[SL][0-9]`, **`.*P[0-9]*$`**, `^C..orf..`, `^[A-Z]{2}[0-9]{5,}`, `\.[0-9]+$`. Note that `.*P[0-9]*$` wrongly removes every normal gene ending in "P + digits" or in "P" (TP53, MMP9, TIMP1, AQP1, FABP4…).
- [code/HF] The released `basecount_1000per_15000max.pkl` is a Python list of **15,012** gene symbols (sorted). Stack-Large and Stack-Large-Aligned ship the same file (925,039 bytes).
- [inference] Implication: genes that are not among the 15,012 (whether because they were never a top-1,000 HVG in any file or because the regex filter removed them, as with TP53, MMP9, TIMP1, IGFBP3, AQP1, BMP4, GBP1) are invisible to the model. This includes, possibly, the very genes targeted by a knockdown; any gene-level evaluation axis should be checked against the vocabulary before using Stack.
- [code] Genes are aligned by exact match on the **upper-cased symbol** (`data/training/datasets.py:723-790`); missing genes are filled with 0; for duplicated names the last one wins (as noted in the tutorial). The authors confirmed in issue #10 that both checkpoints can only predict these 15,012 genes.

### 1.3 Normalisation
- [code] **The input is raw counts.** The model applies `log1p(counts)` internally (`models/core/base.py:153`) **with no CP10k/library-size normalisation** before the log. Library size `l = Σ_g counts` is summed only over the 15,012 model genes (`base.py:151`).
- Output: the decoder gives `softmax(logits) × l` as the NB mean and a `softplus` output as the per-cell, per-gene dispersion θ (`base.py:108-122`).
- The target of the gene-space distribution loss in post-training is `log1p(1e4·x/l)` (log-CP10k, `models/finetune/mixins.py:157-159`).

### 1.4 Cell-set construction and size
- [paper+code] Pre-training: each file is cut, in row order, into **non-overlapping contiguous blocks of K = 256 cells**, and a tail block shorter than K is dropped (`data/training/datasets.py:571-612`). In scBaseCount one SRX is one biological sample; CELLxGENE datasets are mostly sorted by donor.
- [paper] Cell-set size ablation (Fig. S1C, Large, full scBaseCount): among 128/256/512, 256 gives the best validation loss, while reconstruction metrics are best at 128.
- [paper] Fig. 1D: with the total context size fixed and the number of unique cells varied by repetition, Stack only beats the ablation without inter-cell attention once there are >32 unique cells. The larger the model, the more the gain grows with the number of unique cells.
- At inference time (embedding): `TestSamplerDataset` cuts blocks of `sample_size` in file order, and the tail block is filled up by **upsampling with replacement** from the remaining cells (`data/training/datasets.py:864-925`). The authors stress that the input AnnData must be sorted/split by individual biological sample (issue #12).

---

## 2. Tokenization (cell → 100 gene-module tokens)

- [code] `base.py:45-49`: `gene_reduction = Linear(G=15012 → n·d) → GELU → Dropout`, reshaped to `(B, K, n=100, d)` (`base.py:87-90`). Large: d = 16 → 1,600 dims; Base d = 8 (800); XLarge/Huge d = 32 (3,200). In the code `n_hidden` means the token count n = 100. (Section 4.4 of the paper writes "hidden dimension d = 100", which is a typo or a mix-up; the code and the parameter counts in Table 3 show that n = 100.)
- [code] Gene-module identity encoding: `gene_pos_embedding ∈ R^{100×d}` (`base.py:51`, randomly initialised with `randn`) is added to the attention input **in every layer** before intra-cell attention, but **does not enter the residual stream** (`modules/attention.py:112-115`: `cell_attn(x + pos)`, with the residual being `x + attn_out`). Section 4.1.2 of the paper says it is added once to the perceptron output, which differs from the implementation (the code behaves more like an absolute position, without RoPE, injected in every layer).
- There is no external gene-semantics embedding (no ESM, GenePT, etc.); the tokenizer is trained end to end with the model. [paper] In Large, 526 of 699 (75.3%) of the top-10 weighted genes of each token appear in only one token, and GO enrichment shows that tokens are functionally coherent (Fig. 1E/S2).
- No CLS token; the embedding = the 100 tokens of the last layer flattened (1,600 dims).
- Masking: masked genes are simply set to 0 in log1p space (no mask token and no mask-indicator channel) (`base.py:124-140`).

---

## 3. Architecture details

### 3.1 TabularAttentionLayer (`modules/attention.py:61-131`)
Note that the code names are **the opposite** of what they do: `cell_attn` is within a cell (across tokens), and `gene_attn` is across cells.

```
x: (B, K, n, d)
# (1) intra-cell MHA: sequence = 100 tokens, width d
x_cell = x.reshape(B*K, n, d)
x_cell = LN(x_cell + MHA_d(x_cell + pos_emb))        # heads hard-coded to 8 (d=16 → head_dim=2)
# (2) inter-cell MHA: sequence = K cells, width n*d
x_gene = x.reshape(B, K, n*d)
x_gene = LN(x_gene + MHA_{nd}(x_gene, attn_mask))    # heads = n_heads (default 8; paper: XL/Huge = 20)
# (3) token-wise FFN: d → 4d → d, GELU
x = LN(x + FFN(x))
```
- **Post-LN** (LayerNorm after the residual), once per sub-layer.
- Hand-written MHA: `qkv` has no bias, `proj` has a bias; `softmax(QK^T/√h)V`. **Flash attention / `F.scaled_dot_product_attention` is not used**, and there is no varlen, no chunking and no KV cache (`attention.py:28-58`). Boolean masks are filled with `-inf`.
- The intra-cell head count is hard-coded to 8 in the code (`attention.py:76`) and only falls back when it does not divide evenly; the `n_heads` argument only affects inter-cell attention.
- **There is no positional encoding along the cell dimension.** Apart from the query position embedding and the causal mask, the model is therefore equivariant to permutations of cell order, so K can change between pre-training (256) and post-training (512) (`finetune/utils.py:107-115` overrides `n_cells`).
- Dropout: 0.0 in the pre-training config (`configs/training/bc_large.yaml`).
- Decoder: `Linear(nd→2nd) → GELU → Dropout → Linear(2nd → 2G)` (`base.py:67-72`), applied to each cell independently.

### 3.2 Parameter counts (Table 3; reproduced term by term by hand with G = 15,012, exact match)
| Variant | Layers N_L | n·d | Total params | Non-embedding params | Params per layer | Decoder params |
|---|---|---|---|---|---|---|
| Base− | 3 | 800 | 69.1M | 57.0M | 2.6M | 49.4M |
| Base | 6 | 800 | 76.7M | 64.7M | 2.6M | 49.4M |
| Base+ | 9 | 800 | 84.4M | 72.4M | 2.6M | 49.4M |
| Medium | 6 | 1600 | 186M | 163M | 10.2M | 101.2M |
| **Large (released)** | 9 | 1600 | **217M** | **193M** | 10.2M | 101.2M |
| XLarge | 6 | 3200 | 506M | 459M | 41.0M | 212.7M |
| Huge | 9 | 3200 | 629M | 582M | 41.0M | 212.7M |

"Embedding parameters" = the 15,012 × 1,600 of `gene_reduction` ≈ 24.0M. Almost all per-layer parameters are in the inter-cell attention's `qkv (1600→4800)` + `proj (1600→1600)` ≈ 10.24M; intra-cell attention and the FFN together have only ~3k parameters per layer. → **Stack's "depth" is actually shallow: the local computation per token is tiny, and nearly all capacity sits in inter-cell attention and the decoder.**
- Extra parameters in post-training: `query_pos_embedding ∈ R^{100×16}` (1,600 parameters) + the classifier `Linear(3200→100) → GELU → Linear(100→1)` (~320k) (`models/finetune/mixins.py:24-41`).
- The released `.ckpt` files are 2.61 GB (`bc_large.ckpt` 2,610,004,146 B; `bc_large_aligned.ckpt` 2,613,863,242 B), about 3× the fp32 parameter size → [inference] they include the Adam optimiser state.

### 3.3 Masking / batching / efficiency
- Pre-training mask: for each mini-batch a ratio p ~ U(0.1, 0.8) is drawn and ⌊pG⌋ genes are chosen at random and masked **in every sample and every cell of the batch at once** (`base.py:124-140`), i.e. a "rectangular mask" (R²MAE style; the paper's upper bound of 0.8 is above R²MAE's 0.5).
- No attention padding mask: every cell set has a fixed K (tails are filled by upsampling), so varlen is not needed.
- Efficiency comes mainly from data loading: direct h5py reads of CSR, contiguous block reads per file, cached indices. The paper reports ~1.6×10⁴ cells/s, more than 75× faster than State Embedding's loader, and Large pre-training finishes in **2–3 days on a single H100** (bf16 mixed precision, batch 32, 4 workers, 320 GB system RAM). Attention itself is O(K²), but K = 256/512 is small and is not specially optimised.

---

## 4. Training objectives

### 4.1 Pre-training (`base.py:142-204` + `models/core/losses.py`)
- 𝓛 = 𝓛_recon + λ_SW·𝓛_SW, λ_SW = 0.01 (CLI default, `launch_training.py:104`; same in the paper).
- 𝓛_recon: NB negative log-likelihood on the masked genes (`scvi.distributions.NegativeBinomial`, μ = softmax·l, θ from the decoder), averaged over masked entries (`losses.py:14-31`). **Unmasked genes contribute no loss.**
- 𝓛_SW: 32–128 of the K cells are drawn at random (uniform on [32, 128] in the paper), the **set mean is subtracted**, and the sliced Wasserstein distance to a standard-normal sample of the same shape is computed (64 random projections, MSE after sorting) (`losses.py:33-66`, `modules/regularizers.py`). This is equivalent to a prior 𝒩(set mean, I): embedding = a per-cell-set constant offset + a standard normal. The paper uses this to argue for linear identifiability (in the iVAE framework). Ablation (Fig. S1D): removing SW, or removing it together with inter-cell attention, worsens the validation metrics.
- Optimisation: AdamW, peak lr 1e-4 (3e-5 for XL/Huge), wd 3e-3, 1 epoch of linear warmup + cosine, 10 epochs, batch 32 (`configs/training/bc_large.yaml`).

### 4.2 Post-training (ICL alignment) (`models/finetune/model.py`, `mixins.py`, `finetune/lightning.py`, `data/finetuning/datasets.py`)
**Input construction (what the code actually does):**
1. There are two dataset config types (`data/finetuning/datasets.py:35-70`):
   - `human`: group = donor/sample column, identity = cell-type column; query rule `DIFFERENT_GROUP` (same identity, different group, preferring the same file).
   - `drug`: group = condition column, identity = **cell-line** column; query rule `CONTROL_GROUP` (same cell line, condition == control_condition) (`:837-848`, `:927-970`).
2. Each training sample = **K = 512 cells contiguous in row order** within one (file, group) (`:1212-1293`; groups with fewer than K cells are skipped outright).
3. `find_replacement_cells` (`:757-869`): shuffles the identities and picks identities to replace **one whole group at a time** until `replacement_ratio = 0.75` is reached (the last group may be replaced partially). For every replaced cell it finds a query cell of the same identity (a random contiguous interval of the candidate pool, repeated cyclically if it is too short).
4. Class balancing (`:1325-1355`): random identity order; over-represented identities are downsampled to the mean, then everything is upsampled with replacement back to K. Duplicated cells only count towards the distribution loss on their first occurrence (`position_mask`).
5. Returns `ground_truth` (the original cells of the prompt sample, left side) and `observed` (kept cells + query replacement cells, right side).
6. Model forward (`models/finetune/model.py:33-71`):
   - `n_kept_cell = (1−0.75)×512 = 128`: the first 128 positions = **prompt condition**;
   - `time ~ U(0, 0.9375)` is drawn and `n_context = 128 + ⌊384·time⌋`: inputs at positions [128, n_context) are **replaced by the ground-truth cells** (= prompt context / already "unmasked" targets);
   - positions ≥ n_context: query cells (same-identity cells from other samples), with `query_pos_embedding` added;
   - inter-cell attention mask: `attn_mask[:128, 128:] = True` → prompt-condition cells **can only attend to each other**, while all other cells can attend to everything.
   - Gene mask ratio U(0.1, 0.3).
   - **Differences from the paper:** Section 4.2.2 of the paper says the prompt-context cells are "sampled" from the mean/dispersion predicted by the teacher; the code uses the ground-truth counts directly. Also, because cells are sorted by identity after class balancing and the identity order is random, a replaced identity can land in the first 128 "prompt condition" positions (in which case the prompt condition contains query cells). Whether this is intended is **(unverified)**.
7. Teacher (`finetune/lightning.py:109-133`): same architecture; input = the complete original prompt sample (no gene mask, no causal mask, no query embedding); the embeddings at the target positions are the distillation targets. `on_fit_start` copies the weights from the student; **an EMA update with decay 0.95 is applied every 500 steps** (`:43-44, 275-294`); checkpoints only store the student.

**Losses** (`model.py:92-132`, `mixins.py:66-219`):
- 𝓛_FT = 𝓛_recon + 𝓛_dist + λ_SW·𝓛_SW + 𝓛_CLS (λ_recon = λ_CLS = 1; λ_SW takes the checkpoint's sw_weight = 0.01).
- 𝓛_dist = (0.5·𝓛_gene + 0.5·𝓛_embed) **/ (1 − time)** (a masked-diffusion-style weight that the main text of the paper does not mention).
- 𝓛_gene: only on positions ≥ n_context. The top 1,000 variance genes are selected by Pearson residuals with θ = 100 **on the target cells**. Prediction = `ReparamNBLogSampler` (`models/utils.py:118-170`): a reparameterisable approximation of log1p(X/N) under NB(μ, θ) as "zero-inflated log-normal + Concrete gate (τ = 0.7)" (matching the first two moments), with θ set to **the median dispersion of the first n_context positions (detached)**. This is what the paper describes as "using the prompt's dispersion to counter over-smoothing". Target = ground-truth log-CP10k; distance = geomloss `SamplesLoss("energy")`; **stratified by cell type**, counted only for classes with ≥2 cells.
- 𝓛_embed: an energy distance stratified the same way, student embeddings at query positions vs the teacher's target embeddings.
- 𝓛_recon: NB NLL on masked genes, only at positions < n_context (prompt condition + context), with the ground truth as target.
- 𝓛_SW: over all cells, with a fixed subsample of 128.
- 𝓛_CLS: the classifier input is `[mean(prompt-condition embeddings), single-cell embedding]` (all detached); context positions are labelled 0 and query positions 1; BCE with pos_weight balancing. Gradients to the classifier and the query embedding are multiplied by 10 (gradient hook).
- [paper] Post-training hyper-parameters: start from the Large pre-trained weights and train 8 epochs, AdamW lr 2e-5, wd 3e-3, cosine (1 epoch warmup, min 5e-6), batch 8 × grad-accum 4 = 32, a single H100 80GB + 400GB RAM, bf16. Data: the 189 CELLxGENE datasets with >50,000 cells and ≥5 donors (45M cells, with author-annotation columns picked automatically) + Parse 10M PBMC (12 donors, 90 cytokines); train/val/test split by donor/sample. The main text says "~55M".
- [paper] Post-training from scratch (without the pre-trained weights) is clearly worse (Fig. S8; advantage of starting from the pre-trained weights: on Dong cytokines, setting 1, Pearson Δ +11.8%, DE Spearman LFC +20.1%, PR-AUC +15.6%, Overlap@N +35.4%, Jaccard +39.2%, Spearman effect size +5.6%; on OpenProblems +11.5%, +11.7%, +4.5%, +19.2%, +13.4%, **−15.2%**).

---

## 5. In-context inference (code: `models/core/inference.py:551-937`; CLI: `cli/generation.py`)

**Inputs:** `base_adata` (the prompt pool, raw counts), `test_adata` (query cells, raw counts), and a gene list. `--split-column` splits the base by a column (donor / drug / condition), **runs generation once per value**, and writes one h5ad per value.

**Each step** (`get_incontext_prediction`):
1. Base cells are randomly permuted. Each window has `n_cells` positions (= the checkpoint's `n_cells`; Aligned was trained with sample_size = 512 → [inference] 512, not verified by reading the checkpoint hparams), of which `n_test = ⌊n_cells·(1 − prompt − context)⌋` are query cells and the rest are base cells. If there are too few query cells, they are cycled from the start.
2. The first `prompt_ratio·n_cells` base cells = prompt condition (causal mask); the following base cells = prompt context (**real base cells**, not predictions); `query_pos_embedding` is added at query positions. **No gene mask is applied at inference** (`inference.py:159` hard-codes 0).
3. Output: NB mean = softmax × **the query cell's own library size**; `counts = NB(μ, θ = median dispersion over the base positions in the window).sample()` (`:228-236`) → **integer counts**; plus a classifier logit (lower = more prompt-like = higher confidence).
4. Iteration (`get_incontext_generation`; CLI defaults `--mode mdm`, `--num-steps 5`): with T = 5 the mask ratios are [0.8, 0.6, 0.4, 0.2, 0]; the context ratio rises linearly from 0.2 to 0.4; the prompt ratio stays at 0.25. At each step only the fraction of query cells that are **still masked and have the lowest logits** is replaced by predictions, and cells with logit > 0 are then re-masked; at the last step everything is replaced. The replaced counts are written back into test_adata and used as the query input of the next step (**the query input itself is progressively rewritten**).
5. The output is aligned back to test_adata's genes (only 15,012 genes). The final logit is stored in `obs['gen_logit']`; Perturb Sapiens filters out cells with logit > 2.5.
- Known bug: issue #17 reports `np.quantile` failing with "Quantiles must be in the range [0, 1]" (the unmask_rate at `inference.py:700` can go out of range).

**What the prompt is, how many cells, and whether it can be a different perturbation:**
- In training the prompt and the target **come from the same group** (the same donor/sample, or the same condition in the `drug` type), and the query comes from a different group. What the model learns is "transfer the state of the prompt's group onto the query's identity". **What happens when the prompt comes from a "related but different" perturbation is tested neither in the paper nor in the code** → (unverified). [inference] The model would transfer the prompt's own state (i.e. the effect of the related target), not the effect of the intended target.
- Official tutorial (OpenProblems): prompt = **T cells under a given drug**, query = myeloid/B cells under DMSO; each sm_name is generated separately. DMSO itself is also run as a **synthetic control** (the model's prediction for the control condition replaces the real control, to cancel the systematic bias/batch effect introduced by generation).
- Each window's prompt is just ~128 + 102→204 random cells from the base; additional base cells are used by **rotating across windows**, not by putting them into the same context.
- Effect of context composition (embedding task, Fig. S4): restricting the context to the same cell type (CT Context) shrinks the advantage; **randomly shuffling the cell order (Full Shuffle) drops performance to roughly the level of other methods**. Generation task: T = 5 vs T = 1 differs only slightly and in inconsistent directions (Fig. S9; see §6.4).

---

## 6. Experimental results (v1; all numbers from the main text or figure annotations)

### 6.1 Embeddings (zero-shot)
- Linear probes on observational data (Fig. 2D, relative to the best non-Stack method): Kidney disease +62.2% / other +91.7%; Brain (SEA-AD MTG) +8.3% / +0.8%; Lymph node (BCL) +46.4% / +433.5%; LUCA +29.8% / **−2.2%** (LUCA is partly in the training set; State SE does better).
- Linear probes for perturbation classification (Fig. 2E, balanced accuracy, relative to the best alternative): OpenProblems (147 perturbations) +6.0%; Tahoe-100M (250 perturbations from plates 1–3) +139.4%; Parse-PBMC (90) +98.7%; **Xaira / X-Atlas:Orion (top-50 perturbations per cell line, 2 cell lines) +21.2%**.
  - Xaira details (enlarged in Fig. S4): x axis 0–0.04; Stack-Large (BC, All Context) ≈ 0.04 is the highest; second is **PC-HVG**, third State (SE); Stack-Huge/Base, CT-Context and the CxG versions are **all below PC-HVG** (CT Context annotated −15.4%). Chance level is 1/50 = 0.02. The paper itself acknowledges that all methods have very low absolute performance on genetic-perturbation classification, possibly because of measurement noise and because the expression profiles of the perturbations are similar to one another.
- scIB integration: first on all 4 observational datasets, +1.8% over the best State (SE); first on 21 of 25 Tabula Sapiens tissues.

### 6.2 Cell prompting / ICL (Fig. 3; post-trained Large, T = 5; percentages = mean improvement over the best non-Stack baseline)
| Setting | Data | Pearson Δ | DE Spearman LFC | PR-AUC | DE Overlap | Other |
|---|---|---|---|---|---|---|
| 1 perturbation → new cell type (same sample) | Dong 2023 cytokines (6) | +45.8% | +3.5% | +55.5% | +5.2% | Spearman effect size +23.0% |
| 1 | OpenProblems drugs (12) | +11.1% | +3.5% | +10.0% | +119.8% | Spearman effect size +8.9% |
| 2 perturbation → new sample (cross-dataset, T cells) | Parse prompt → Dong control (7 conditions) | **−0.7%** | **−14.2%** | +37.7% | +25.4% | Precision@N +45.2% |
| 3 observational: held-out cell type (same dataset, across donors) | 4 atlases | DE direction +1.5% | **−21.7%** | +14.7% | +13.7% | Spearman effect size +27.2% |
| 4 cross-dataset cell-type generation | 4 PBMC prompts + Tabula Sapiens query | DE direction +22.8% | +25.7% | +4.6% | +5.1% | normalised scIB +94.5%, Spearman effect size +41.6% |
- First on 28 of 31 evaluations in total. Baselines: the raw query; the nearest/same cell type in the prompt; PerturbMean/DonorMean (an oracle in settings 3/4); scVI; State (ST+SE, **trained specifically** on prompt + query + auxiliary data; latent dim 328, set length 32, 60k steps).
- Metrics (cell-eval v0.6.6, Wilcoxon + BH, top-2000 log-normalised HVGs): Pearson Δ, DE Spearman LFC, PR-AUC, DE overlap accuracy, Spearman effect size, DE direction match, Precision@N, Jaccard.
- Failures/weaknesses (as stated by the paper): in setting 2, because Parse and Dong differ in stimulation dose/duration, **cross-dataset transfer of weak-effect cytokines such as TNF-α and IL-6 is close to zero for all methods**; "nearest/same cell type in the prompt" is a strong baseline; TNF-α in epithelium goes in the **opposite** direction to the in-vitro data (Table S1: Pearson Δ = −0.122***), which the authors attribute to secondary signalling. The Discussion explicitly lists calibration for rare cell types and **weak perturbation effects** as unsolved, and notes that post-training is dominated by in-vivo immune cells, which limits generalisation to in-vitro perturbation studies.

### 6.3 Perturb Sapiens (v1)
- Prompt: 90 cytokines from Parse donor 1 + 111 drugs from OpenProblems donor 2 (the Fig. 4A caption says 144 drugs; the main text says 111). Query: tissue-balanced Tabula Sapiens. 513,870 cells per condition, 28 tissues, 40 broad cell classes, 201 perturbations. The HF dataset is 1.34 TB, CC BY 4.0. The v2 abstract says it was extended to 892 perturbations (including genetic ones) — **(unverified)**.

### 6.4 Ablations
- Generation steps T = 5 vs T = 1 (Fig. S9, T = 5 relative to T = 1): Dong (setting 1) Pearson Δ −8.0%, DE Spearman LFC −14.9%, PR-AUC +3.5%, DE overlap +17.3%, Jaccard +13.6%, Spearman effect size −1.7%; OpenProblems −2.1%, +4.5%, +2.5%, +6.3%, −2.6%, +12.5%; setting 2 Pearson Δ +12.5%, DE Spearman LFC +16.7%, PR-AUC +0.8%, DE overlap −8.2%, Jaccard −7.0%, Precision@N −8.8%; settings 3/4 all change by less than ±6% (setting 3: Spearman effect size −1.4%, Pearson Δ +5.9%). → The benefit of multi-step generation is **small and inconsistent in direction**.
- Scale (Fig. S1A/B): a larger latent dimension gives steady gains; going deeper gives similar validation loss on full scBaseCount but improves other metrics.
- Pre-trained weights vs post-training from scratch: see §4.2 (Fig. S8).
- There is **no** ablation of prompt size or prompt composition (e.g. mixed conditions, related conditions) on ICL generation (v1).

---

## 7. Code and engineering

### 7.1 Repository layout (`src/stack/`, ~9.3k lines of Python)
- `modules/attention.py` (tabular block), `modules/regularizers.py` (SW)
- `models/core/base.py` (tokenizer/layers/decoder/pre-training forward), `losses.py`, `inference.py` (embedding, prediction, ICL generation)
- `models/finetune/model.py`, `mixins.py` (post-training forward and losses); `models/utils.py` (ReparamNBLogSampler, gene alignment)
- `data/training/datasets.py` (pre-training / inference datasets), `data/finetuning/datasets.py` (post-training episode sampling, 2,283 lines), `data/hvg.py`, `data/gene_processing.py`, `data/h5_manager.py`
- `training/`, `finetune/` (Lightning modules and datamodules), `cli/` (4 entry points)
- `configs/training/bc_large.yaml`, `configs/finetuning/ft_parsecg.yaml`; `notebooks/tutorial-embed.ipynb`, `tutorial-predict.ipynb`; `tests/`
- Package name `arc-stack` 0.1.3 (PyPI). Dependencies: anndata, h5py, **geomloss**, numpy, pandas, psutil, pytorch-lightning ≥ 2.1, PyYAML, scipy, **scvi-tools** (only for the NB distribution), torch ≥ 2.0, wandb. README test environment: Ubuntu 22.04, Python 3.10.18, PyTorch 2.5.1+cu121, H100 80GB.
- Console entry points: `stack-train`, `stack-finetune`, `stack-embedding`, `stack-generation`.

### 7.2 Released checkpoints
| HF repo | File | Size | Notes |
|---|---|---|---|
| `arcinstitute/Stack-Large` | `bc_large.ckpt` | 2.61 GB | Pre-trained, 217M, full scBaseCount, 10 epochs |
| `arcinstitute/Stack-Large-Aligned` | `bc_large_aligned.ckpt` | 2.61 GB | + CELLxGENE 45M + Parse 10M, 8 epochs |
| Both repos | `basecount_1000per_15000max.pkl` | 925 KB | 15,012 genes |
- Not gated; can be downloaded directly. Only Large is released (Base/Medium/XL/Huge are not). There are also example datasets `arcinstitute/Stack-DrugPBMC-Example` and `arcinstitute/Stack-CellxGene45M`.

### 7.3 Compute requirements
- Paper: both pre-training and post-training ran on a single H100 80GB (with 320GB and 400GB of system RAM respectively); pre-training takes 2–3 days; post-training time is **not reported (unverified)**.
- Inference: 217M parameters, K = 512, batches of 32 windows → [inference] a single 24–40 GB GPU should be enough (not measured). Nothing was executed for this review.

### 7.4 Licence
- Code: `LICENSE` = **CC BY-NC-SA 4.0** (note that `setup.cfg` says `license = Apache-2.0`, which contradicts the LICENSE file; LICENSE/README should be taken as authoritative).
- Weights and outputs: Arc Research Institute Stack Model Non-Commercial License (2025-06-23) + Acceptable Use Policy. Key points: only Non-Commercial Purposes are allowed; research initiated or funded by a Commercial Entity does **not** count as non-commercial; Derivative Work includes fine-tuning, adjusting weights, distillation/pruning, continued training, **and "use as a prior"**; Output (data generated by the model) is also covered, and redistribution must include the licence and a citation; commercial use requires contacting legal@arcinstitute.org. → No commercial use; commercial or company-affiliated use requires confirmation (**unverified**).
- This summary is a technical reading of the licence texts and is not legal advice. Consult the licence documents themselves (and, if needed, a lawyer) before relying on it.

### 7.5 Running in-context prediction on your own data (no training)
```bash
pip install arc-stack   # needs a GPU machine
huggingface-cli download arcinstitute/Stack-Large-Aligned --local-dir Stack-Large-Aligned
# base.h5ad: prompt cells (raw counts; var index, or --gene-name-col, holds HGNC symbols), obs['pert'] = condition
# test.h5ad: a set of control cells from the new context (raw counts); one prediction is produced per control cell
stack-generation --checkpoint Stack-Large-Aligned/bc_large_aligned.ckpt \
  --base-adata base.h5ad --test-adata test.h5ad \
  --genelist Stack-Large-Aligned/basecount_1000per_15000max.pkl \
  --split-column pert --split-values TARGET_X non-targeting \
  --prompt-ratio 0.25 --context-ratio 0.4 --num-steps 5 --batch-size 16 \
  --output-dir out/
```
- Each split value yields one h5ad: every query control cell maps to one "predicted post-perturbation" counts cell (15,012 genes, obs aligned to test).
- Using this to predict perturbed cells for a set of control cells of a new context [inference]: (i) optionally run `non-targeting` as the prompt to get a synthetic control, compute Δ = pred(X) − pred(ctrl), and add it to the **real** controls in log/count space, to cancel the generator's systematic bias; (ii) genes outside the 15,012-gene vocabulary are not predicted and have to be carried over from the controls unchanged; (iii) if a knocked-down target is not in the vocabulary, its knockdown has to be written in by a separate rule; (iv) library size comes from the query cell's total over the 15k genes, so output depth follows the control cells.

### 7.6 Post-training on custom episodes
- Ready-made route (cross-cell-line transfer of the same target): put the CRISPRi data into an h5ad (raw counts; `obs` has a target column and a cell_line column; all cell lines in **one file**, so that intra-file replacement can find same-cell-line controls), then:
```bash
stack-finetune --checkpoint_path Stack-Large/bc_large.ckpt \
  --dataset_configs "drug:/data/crispri_dir:target_gene:cell_line:non-targeting:false:gene_symbol" \
  --genelist_path basecount_1000per_15000max.pkl --sample_size 256 --replacement_ratio 0.75 \
  --batch_size 8 --accumulate_grad_batches 4 --learning_rate 2e-5 --max_epochs 8
```
  Each episode: 256 contiguous cells from one target group → some cell lines are chosen at random and wholly replaced by non-targeting cells of the same cell line → the model learns the knockdown effect in those cell lines from the same-target cells of the cell lines that were not replaced.
- Pitfalls (from reading the code):
  1. Any (file, target) with fewer than `sample_size` cells is skipped entirely (`datasets.py:1252`). Typical CRISPRi screens have only a few hundred cells per target per cell line, so use 128/256.
  2. When a group contains only one cell line, cells of that same cell line are partially replaced, which means the prompt contains perturbed cells of the same cell line → **leakage**; this does not represent zero-shot cross-cell-line prediction. `find_replacement_cells` has to be rewritten so that perturbed cells of the target cell line never appear in the prompt.
  3. The train/val/test split is by group (= target) (`_split_groups`), not by cell line; held-out evaluation has to be by cell line.
  4. "Pathway-related targets as the prompt" needs modification: redefine the group as "target X + a set of related targets", build prompts with mixed conditions, and **add perturbation-identity conditioning** (e.g. add an embedding of the target gene at query positions, similar to `query_pos_embedding`); otherwise the model has no way of knowing that it should predict X.
  5. The gene mask defaults to U(0.1, 0.3), and the distribution loss uses only the top 1,000 HVGs. For weak-effect targets, the target gene itself and small-effect DEGs may not be among those 1,000, so either force-include the target genes or switch to all genes.

---

## 8. Assessment for cross-context perturbation prediction

### 8.1 What Stack already does (overlap with the approach under consideration)
- Combines control cells and other context cells into one set and performs in-context prediction with **bidirectional inter-cell attention**; causal mask on the prompt block, query position embedding, iterative masked-diffusion generation, confidence classifier.
- The code already has episode sampling with "same-cell-line controls as the query, the same condition in other cell lines as the prompt" (the `drug` type).
- Outputs counts directly (NB sampling) and preserves the query cells' library size.

### 8.2 What Stack lacks
1. **Genetic perturbations**: v1 post-training has no genetic-perturbation data; the only genetic result is the Xaira probe (absolute values near chance, conclusions not robust); there is no ICL validation on weak-effect CRISPRi. v2 claims to add genetic perturbations — **(unverified)**.
2. **Pathway-related prompts**: both training and inference assume prompt = the same condition; there is no mechanism for "related but different perturbations", and no perturbation-identity input or gene-knowledge prior.
3. **Zero-shot across cell lines**: no experiment holds out cell lines; the main evaluations are on PBMC/in-vivo immune cells.
4. **Gene axis**: a fixed 15,012-gene HVG union passed through a gene-module projection; genes outside it (including regex casualties such as TP53, and possibly knocked-down targets) are invisible; the vocabulary is selected by HVG and may favour highly variable genes over weak-effect DEGs.
5. **Count output**: available, but on 15k genes, NB-sampled with the prompt's median dispersion, and not calibrated to any particular evaluation metric; the output carries the generator's systematic bias (the authors themselves need a synthetic control).
6. Efficiency: hand-written attention, no flash attention, fixed K. Not a big issue at Stack's scale, but going to a larger K (e.g. a 1–2k-cell context) would call for SDPA.

### 8.3 Concrete integration routes (from lowest to highest effort)
- **A. Frozen Stack as an encoder (no training):** feed the new cell line's controls (and the prompt cells) through Stack-Large to extract 1,600-dim context-aware embeddings, and use them as conditioning/features for a flow-matching model. Risks: licence ("use as a prior" counts as a Derivative Work); 15,012-gene vocabulary coverage; weak discriminative power of the embeddings for genetic perturbations (Xaira result).
- **B. Stack-Aligned zero-shot + synthetic-control Δ as a baseline/ensemble member:** only meaningful for targets that appear in other cell lines in the training set; expected to do poorly on weak-effect CRISPRi (setting 2 in the paper already shows that weak effects transfer at close to zero across protocols). A cheap comparison point, but low priority.
- **C. Continue post-training Stack-Large on CRISPRi episodes:** reuse the `drug`-type sampler, but rewrite the replacement logic (avoid leakage, hold out by cell line), add target conditioning, add target genes to the loss gene set, and handle out-of-vocabulary targets (adding them to the vocabulary changes the dimensions of `gene_reduction`/the decoder, which means re-initialising these two layers, ≈125M parameters). Feasible on a single H100, but both the workload and the licence risk are high.
- **D. Re-implement the tabular block inside a flow-matching model (recommended):** ~130 lines are enough to replicate it (intra-cell token MHA + inter-cell MHA over the flattened n·d + token FFN, post-LN), using `F.scaled_dot_product_attention` with a custom block mask: new-cell-line control block / related-target prompt block / noisy query block. Use Stack's "causal prompt block + query position embedding + stratified energy distance + time weighting" as design references, and work directly on the project's own gene axis with a custom tokenizer (e.g. gene-module projection + an explicit target-gene token). No dependence on Arc's weights at all → the licence issue is avoided. Engineering details worth borrowing from Stack: rectangular gene mask, SW regulariser (set-centred + N(0, I)), prompt median dispersion, class-balanced sampling, contiguous-block h5 reads.

### 8.4 Risk register
| Risk | Description | Severity |
|---|---|---|
| Licence | Code CC BY-NC-SA; weights/outputs non-commercial, with Derivative Work covering fine-tuning and "use as a prior"; commercial or company-affiliated use requires confirmation | High (needs legal confirmation, unverified; not legal advice) |
| Weak genetic-perturbation ability | No genetic post-training in v1; Xaira probe ~0.04 vs 0.02 chance, only 1 Stack variant beats PC-HVG | High |
| Gene vocabulary mismatch | Fixed 15,012-gene vocabulary; out-of-vocabulary genes (possibly including knocked-down targets) are invisible; regex wrongly removes TP53 and others | High |
| Weak effects | The paper itself says weak-effect calibration is unsolved; cross-protocol transfer of weak effects is near zero | High |
| Prompt semantics | Only trained with "same-condition prompts"; behaviour with pathway-related prompts is unknown | Medium–high |
| Generation bias | Needs a synthetic control; multi-step generation gains are unstable | Medium |
| Compute | Inference is cheap; post-training needs a single H100 + ~400GB RAM; 2.6 GB checkpoint | Low–medium |
| Code quality | Discrepancies between paper and code (how prompt context is sampled, how the position embedding is injected, K_kept naming); quantile bug in generation (#17); wrong licence field in `setup.cfg` | Low–medium |
| v2 not read | v2 (2026-06) may already add genetic perturbations and change the method | Medium (recommend obtaining the v2 PDF manually and re-checking) |

---

## Sources
- bioRxiv v1/v2 pages: https://www.biorxiv.org/content/10.64898/2026.01.09.698608v1 (v2: …v2; full text blocked, only metadata/abstract obtained via api.biorxiv.org: https://api.biorxiv.org/details/biorxiv/10.64898/2026.01.09.698608 )
- PMC full text (v1): https://pmc.ncbi.nlm.nih.gov/articles/PMC12803207/
- Supplementary PDF + main figures: Europe PMC supplementary-files API https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12803207/supplementaryFiles
- Code: https://github.com/ArcInstitute/stack @ cacc2e4; GitHub issues #10, #12, #14, #15, #17, #18
- Weights: https://huggingface.co/arcinstitute/Stack-Large , https://huggingface.co/arcinstitute/Stack-Large-Aligned (README, file tree, gene-list pkl)
- Data: https://huggingface.co/datasets/arcinstitute/Perturb-Sapiens
- Tool page / press release: https://arcinstitute.org/tools/stack , https://arcinstitute.org/news/foundation-model-stack
