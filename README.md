# A Multi-Omics Integrative Framework Reveals a 9-Gene Prognostic Signature in Triple-Negative Breast Cancer

> **Note:** This repository hosts the analysis code and resources for Abhishek S R’s MSc Bioinformatics dissertation project at Garden City University. It is designed to be a clean, reproducible, and extensible framework for TNBC prognostic biomarker discovery using multi-omics integrative analysis.

---

## 1. Project overview

Triple-negative breast cancer (TNBC) is an aggressive, heterogeneous breast cancer subtype defined by the absence of estrogen receptor (ER), progesterone receptor (PR), and HER2 expression, accounting for ~10–20% of all breast cancers yet contributing disproportionately to mortality.[file:1]
Because TNBC lacks approved targeted therapies for the entire patient population, there is a critical need for robust prognostic biomarkers that can stratify patients into clinically meaningful risk groups.[file:1]

This project develops an integrated *multi-omics* and *network-driven* framework to identify and validate a **9-gene prognostic signature** for TNBC using RNA-seq data from TCGA and independent GEO cohorts.[file:1]
The pipeline combines differential expression, co-expression network analysis, protein–protein interaction (PPI) centrality, survival modeling (LASSO-Cox), and immune microenvironment deconvolution.[file:1]

### Key objectives

- Identify differentially expressed genes (DEGs) between TNBC tumors and adjacent normal tissues.
- Integrate DEGs with co-expression modules to obtain a high-confidence gene set.
- Construct a PPI network and prioritize hub genes using network centrality.
- Build a robust LASSO-Cox prognostic model and derive a stable multi-gene risk signature.
- Validate the signature across independent GEO cohorts without model re-fitting.
- Characterize associations between the signature and the tumor immune microenvironment using CIBERSORT.[file:1]

---

## 2. Data and resources

### 2.1 Discovery cohort (TCGA-BRCA TNBC subset)

- **Source:** TCGA Breast Invasive Carcinoma (BRCA) RNA-seq dataset.[file:1]
- **Curation:** TNBC samples curated by cross-referencing clinical metadata from cBioPortal, retaining only ER-, PR-, and HER2-negative cases.[file:1]
- **Final cohort:** 115 TNBC tumor samples and 113 adjacent solid normal tissue samples.[file:1]
- **Pre-processing:**
  - Batch effect correction using ComBat (sva R package).
  - Quality control using PCA and sample-wise correlation matrices.[file:1]

### 2.2 External validation cohorts (GEO)

Two independent RNA-seq GEO datasets are used purely for *external validation* of the prognostic signature:

- **GSE142258** (n = 15)
- **GSE142731** (n = 42)

These cohorts are never used in model training; they assess generalizability and cross-cohort robustness.[file:1]

> **Dataset schematic / flow diagram**  
> ![Dataset overview](results/plots/risk_score_distribution.png)

---

## 3. Methods and analysis pipeline

The project is implemented primarily in **R**, with standard Bioconductor and CRAN packages for RNA-seq analysis, network construction, survival modeling, and immune deconvolution.[file:1]
Below is a conceptual breakdown of each step; you can mirror this in the repository structure (e.g., under `scripts/`).

### 3.1 Differential gene expression analysis (DESeq2)

- **Goal:** Identify TNBC vs normal DEGs at the gene level.
- **Tool:** DESeq2 (v1.38.0).
- **Filter:**
  - Protein-coding genes only.
  - Significant DEGs defined as:  
    - |log2 fold change| ≥ 1  
    - Adjusted p-value (FDR) ≤ 0.05.[file:1]
- **Outputs (suggested files):**
  - `results/deseq2_differential_expression.tsv`
  - `results/plots/Volcano_plot_DESeq2.png`
  - `results/plots/Top50_DEG_heatmap_DESeq2.png`

> **DESeq2 plots**  
> - Volcano plot: ![DESeq2 volcano](results/plots/Volcano_plot_DESeq2.png)  
> - MA and dispersion: ![MA plot](results/plots/MA_plot_DESeq2.png), ![Dispersion plot](results/plots/Dispersion_plot_DESeq2.png)  
> - Top 50 DEG heatmap: ![Top 50 DEGs](results/plots/Top50_DEG_heatmap_DESeq2.png)

### 3.2 Functional enrichment (GO / KEGG)

- **Tool:** clusterProfiler (v4.6.0).
- **Analyses:**
  - GO Biological Process (BP), Molecular Function (MF), Cellular Component (CC).
  - KEGG pathway enrichment.
- **Highlights:** DEGs enriched in processes such as small GTPase-mediated signal transduction, cell adhesion, nucleotide metabolism, and cell cycle-related pathways.[file:1]
- **Outputs (suggested):**
  - `results/go_enrichment_degs.tsv`
  - `results/kegg_enrichment_degs.tsv`
  - `results/plots/GO_BP_dotplot_publication.png`
  - `results/plots/KEGG_dotplot.png`

> **Enrichment plots**  
> - GO BP dotplot: ![GO BP](results/plots/GO_BP_dotplot_publication.png)  
> - KEGG dotplot: ![KEGG](results/plots/KEGG_dotplot.png)

### 3.3 Co-expression network analysis (CEMiTool)

- **Tool:** CEMiTool (v1.24.0).
- **Input:** VST-normalized TNBC expression matrix.
- **Steps:**
  - Automatic soft-threshold selection, enforcing approximate scale-free topology.[file:1]
  - Identification of gene co-expression modules.
  - Gene set enrichment analysis (fgsea) to identify modules associated with tumor vs normal groups.[file:1]
- **Key result:** Module **M3** shows the strongest positive association with the TNBC tumor group and is enriched in cell cycle and mitotic signaling pathways.[file:1]
- **Intersection for downstream analysis:**
  - DEGs ∩ M3 module genes → 1,030 high-confidence genes.

> **Module / intersection figure**  
> *(You can export and save your CEMiTool module and Venn plots under `results/plots/` and embed here similarly.)*

### 3.4 PPI network construction and hub gene prioritization

- **Tool:** STRING for PPI, CINNA for centrality analysis.[file:1]
- **Steps:**
  - Build PPI network from the 1,030 intersected genes.
  - Extract giant component: 581 nodes, 8,060 edges.[file:1]
  - Apply network centrality measures via CINNA with PCA to identify the most informative centrality metrics (Closeness Latora and Harmonic).[file:1]
  - Select 266 hub genes with centrality scores above the mean.[file:1]
- **Examples of top hub genes:** IL6, ESR1, FN1, BCL2, MMP9, IGF1, PPARG, EGF, CXCL12.[file:1]

> **PPI / centrality plots**  
> *(Embed your STRING/CINNA plots here once exported, e.g.: `results/plots/ppi_network.png`, `results/plots/hub_centrality.png`.)*

### 3.5 Survival analysis and LASSO-Cox model

- **Goal:** Derive a robust prognostic gene signature from hub genes.
- **Tools:**
  - `survival` R package (Cox regression).
  - `glmnet` (v4.1.6) for LASSO-penalized Cox modeling.[file:1]
- **Workflow:**
  1. Univariate Cox regression on 266 hub genes.  
     - Genes with p ≤ 0.10 carried forward (28 candidates; 15 at p ≤ 0.05).[file:1]
  2. LASSO-Cox with 10-fold cross-validation at the optimal λ (λ_min = 0.019).[file:1]
  3. Consensus LASSO over 100 bootstrap iterations, retaining only genes selected in ≥50% of iterations to stabilize against low EPV (events-per-variable). [file:1]

- **Final 9-gene prognostic signature:**
  - **Genes:** IL12RB2, DMD, RAD51, KLF4, OXTR, BMP7, VWF, SERPINE1, COL9A3.[file:1]
  - Negative coefficients (protective): IL12RB2, DMD, OXTR, COL9A3, BMP7.[file:1]
  - Positive coefficients (risk): RAD51, KLF4, SERPINE1, VWF.[file:1]

- **Performance in TCGA TNBC cohort:**
  - Apparent C-index ≈ 0.899.
  - Optimism-corrected C-index ≈ 0.851 (95% CI ≈ 0.786–0.905).[file:1]
  - Risk stratification based on median risk score: clear separation of high- vs low-risk groups with highly significant log-rank p-value (≈ 5.38 × 10⁻⁸).[file:1]

> **LASSO and survival plots**  
> - Cross-validation curve: ![CV curve](results/plots/cv_curve_filtered.png)  
> - Lambda stability: ![Lambda stability](results/plots/lambda_stability.png)  
> - Selection frequency: ![Selection frequency](results/plots/selection_frequency_barplot.png)  
> - Final coefficients: ![Final coefficients](results/plots/final_coefficients.png)  
> - KM curve: ![Kaplan–Meier](results/plots/kaplan_meier_curve.png)  
> - Time-dependent ROC / AUC: ![Time-dependent ROC](results/plots/time_dependent_ROC_curves.png), ![AUC over time](results/plots/AUC_over_time.png)  
> - Calibration / bootstrap: ![Calibration](results/plots/calibration_plot.png), ![Bootstrap](results/plots/bootstrap_validation.png)  
> - Individual ROC curves: ![Individual ROC](results/plots/individual_ROC_curves.png)

### 3.6 External validation of the 9-gene signature

- **Cohorts:** GSE142258 (n = 15) and GSE142731 (n = 42), merged into a 57-sample validation set.[file:1]
- **Pre-processing:**
  - Restrict to shared genes across cohorts.
  - Minimum count filter: ≥10 counts in ≥10% of samples.
  - Variance-stabilizing transformation (VST) via DESeq2 with `blind = TRUE`.[file:1]
- **Validation strategy:**
  - Apply fixed LASSO-Cox coefficients from TCGA model (no re-fitting).
  - Use the *training-derived* median risk cutoff (median ≈ 4.773) to classify high- vs low-risk.[file:1]

- **Key validation findings:**
  - Risk scores range from ≈ −6.72 to 0.51 (mean ≈ −3.81 ± 1.18).[file:1]
  - Approximately 79% patients classified as high-risk, 21% as low-risk.[file:1]
  - Bimodal and monotonic risk score distribution, consistent with the training model behavior.[file:1]
  - Clear separation of risk groups within each cohort (stratified boxplots, no overlap).
  - Cross-cohort Wilcoxon p ≈ 0.59, indicating no cohort-specific bias in risk distribution.[file:1]

> **Validation plots**  
> - Risk score distribution: ![Risk distribution](results/plots/risk_score_distribution.png)  
> - Risk vs outcome: ![Risk vs outcome](results/plots/risk_distribution_by_outcome.png)

### 3.7 Immune microenvironment analysis (CIBERSORT)

- **Tool:** CIBERSORT (via IOBR R package v0.99.8).[file:1]
- **Input:** TCGA TNBC expression data.
- **Goal:** Estimate relative abundances of 22 immune cell subtypes and relate them to risk groups and signature genes.[file:1]
- **Analyses:**
  - Stacked barplots and heatmaps of immune cell fractions across Normal, Tumor, Low-risk, and High-risk groups.[file:1]
  - Violin plots for differential infiltration patterns.
  - Spearman correlations between 9-gene signature expression and immune cell fractions.[file:1]

> **Immune infiltration plots**  
> - Violin / boxplots: ![Immune violin](results/plots/01_Violin_Boxplot.png)  
> - Stacked barplots (per sample, per group): ![Per sample](results/plots/02_Stacked_Bar_per_Sample.png), ![Mean per group](results/plots/03_Mean_Stacked_Bar_by_Group.png)  
> - Boxplot grid by cell type: ![Boxplot grid](results/plots/04_Boxplot_Grid_by_CellType.png)  
> - Spearman correlation heatmap: ![Immune correlation](results/plots/05_Spearman_Correlation_Heatmap.png)  
> - Sample-wise immune heatmap: ![Immune heatmap](results/plots/06_Sample_Immune_Heatmap.png)  
> - Additional immune plots: dotplot, ridge plot, Wilcoxon tests, etc.:  
>   - ![Immune dotplot](results/plots/07_Dot_Plot_Immune.png)  
>   - ![Immune ridge plot](results/plots/08_Ridge_Plot_Immune.png)  
>   - ![Tumor vs normal](results/plots/09_Tumor_vs_Normal_Wilcoxon.png)  
>   - ![Low vs high risk](results/plots/10_Low_vs_High_Risk_Wilcoxon.png)  
>   - ![CIBERSORT p-values](results/plots/11_CIBERSORT_Pvalue_Distribution.png)  
>   - ![Lollipop TN vs N](results/plots/12_Lollipop_Tumor_vs_Normal.png)  
>   - ![Lollipop low vs high](results/plots/14_Lollipop_Low_vs_High_Risk.png)

---

## 4. Repository structure (suggested)

Below is a suggested layout for this repository. You can adjust as needed.

```text
tnbc-multiomics-prognostic-signature/
├── data/
│   ├── tcga_brca_tnbc/           # Metadata and processed expression (no raw TCGA data in GitHub)
│   ├── geo_gse142258/
│   └── geo_gse142731/
├── scripts/
│   ├── 01_data_preprocessing.R
│   ├── 02_deseq2_differential_expression.R
│   ├── 03_go_kegg_enrichment.R
│   ├── 04_cemitool_coexpression.R
│   ├── 05_string_ppi_cinna_centrality.R
│   ├── 06_lasso_cox_prognostic_model.R
│   ├── 07_external_validation_geo.R
│   └── 08_cibersort_immune_infiltration.R
├── results/
│   ├── tables/
│   └── plots/
├── images/                      # Final figures for README and manuscript
├── docs/
│   └── methods_notes.md         # Extra methodological details if needed
├── poster/
│   └── TNBC_9gene_poster.pdf    # Optional: conference / defense poster
├── LICENSE
└── README.md
```

> **Important:** Do **not** upload raw TCGA or restricted GEO data directly to GitHub. Instead, provide scripts and instructions to download and pre-process data from the original sources.

---

## 5. Reproducibility and environment

You can document and automate the environment using either `renv` or a Docker image.

### 5.1 R package environment (example)

List key dependencies in a file such as `docs/packages_used.md` or manage them via `renv`:

- DESeq2
- sva
- clusterProfiler
- CEMiTool
- CINNA
- survival
- glmnet
- IOBR
- corrplot
- ggplot2
- genekitr

### 5.2 Example usage pattern

A typical analysis run might follow this sequence:

```bash
# 1. Preprocess TCGA and GEO data
Rscript scripts/01_data_preprocessing.R

# 2. Run DESeq2 and export DEG results and plots
Rscript scripts/02_deseq2_differential_expression.R

# 3. Run functional enrichment
Rscript scripts/03_go_kegg_enrichment.R

# 4. Build co-expression modules with CEMiTool
Rscript scripts/04_cemitool_coexpression.R

# 5. Construct PPI and run centrality analysis
Rscript scripts/05_string_ppi_cinna_centrality.R

# 6. Train LASSO-Cox model and derive 9-gene signature
Rscript scripts/06_lasso_cox_prognostic_model.R

# 7. Apply model to GEO cohorts for external validation
Rscript scripts/07_external_validation_geo.R

# 8. Run CIBERSORT-based immune infiltration analysis
Rscript scripts/08_cibersort_immune_infiltration.R
```

You can refine these scripts into reproducible workflows (e.g., `targets`, `drake`, Snakemake`, or Nextflow) as a future extension.

---

## 6. Biological interpretation and significance

- The 9-gene signature captures both intrinsic tumor biology (cell cycle, DNA repair, extracellular matrix remodeling) and the state of the tumor immune microenvironment (NK cells, macrophage polarization, T-cell activity).[file:1]
- High-risk TNBC patients exhibit transcriptional and immune features consistent with more aggressive disease and poorer prognosis.[file:1]
- The optimism-corrected C-index and clean validation across independent cohorts support the potential clinical utility of this signature for risk stratification in TNBC.[file:1]

You can expand this section with more detailed gene-level annotations and literature links in a separate `docs/biological_interpretation.md` file.

---

## 7. Limitations and future directions

Summarizing from the dissertation and providing context for further work:[file:1]

- The model is trained primarily on bulk RNA-seq; it does not account for intratumoral heterogeneity at the single-cell level.
- Prospective clinical validation and integration with clinical covariates (e.g., stage, grade, treatment regimen) are required before translation into practice.
- Functional validation of individual genes (e.g., SERPINE1, RAD51, IL12RB2) in experimental TNBC models is essential to confirm causal roles.
- Future extensions could integrate other omics layers (e.g., somatic variants, methylation, proteomics) and leverage more advanced machine learning frameworks.

---

## 8. How to cite

If you use this repository or the 9-gene TNBC prognostic signature in your work, please cite:

> Abhishek S R, *A Multi-Omics Integrative Framework Reveals a 9-Gene Prognostic Signature in Triple-Negative Breast Cancer*, MSc Bioinformatics Dissertation, Garden City University, 2026.

Once your dissertation is deposited or a manuscript is published, you can update this section with DOI / journal information.

---

## 9. Contact

- **Author:** Abhishek S R  
- **Affiliation:** MSc Bioinformatics, Department of Life Sciences, Garden City University, Bengaluru, India  
- **GitHub:** https://github.com/Abhishek-MSBI  
- **Email:** *(add your preferred contact email here)*

Feel free to open issues or pull requests if you adapt this framework, extend the analysis, or identify any bugs.
