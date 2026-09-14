# Integrative Identification and Independent Validation of miRNA Biomarkers Associated with Triple-Negative Breast Cancer

## Project Overview

Triple-negative breast cancer (TNBC) is an aggressive breast cancer subtype characterized by the absence of estrogen receptor (ER), progesterone receptor (PR), and HER2 expression. This project uses an integrative bioinformatics workflow to identify and independently validate candidate microRNAs (miRNAs) associated with TNBC using publicly available TCGA-BRCA and GEO datasets.

The analysis progresses from discovery and candidate selection to independent validation, ROC analysis, survival analysis, and target-gene investigation.

## Research Question

Can differentially expressed miRNAs be identified as candidate biomarkers associated with triple-negative breast cancer and independently validated across multiple publicly available datasets?

## Objectives

1. Identify differentially expressed miRNAs in clinically defined TNBC.
2. Compare TNBC with normal breast tissue, HR+/HER2− breast cancer, and HER2+ breast cancer.
3. Prioritize candidates based on statistical significance, effect size, and consistency of expression.
4. Independently validate candidate miRNAs using GEO datasets.
5. Evaluate the TNBC discriminatory potential of the leading candidate using ROC/AUC analysis.
6. Assess the association between candidate miRNA expression and overall survival.
7. Investigate experimentally supported target genes and their functional context.

## Study Workflow

TCGA-BRCA miRNA Expression → Clinical Annotation (ER/PR/HER2) → Clinical Subgroup Classification → Differential Expression (DESeq2) → Candidate Selection (padj < 0.05, |log2FC| ≥ 1, direction consistency) → Independent GEO Validation (GSE40049 and GSE86277) → miR-934 Prioritization → ROC/AUC + Survival + Target Analysis → Candidate TNBC Classification Biomarker

## Datasets

### TCGA-BRCA

TCGA-BRCA miRNA expression data were obtained through the Genomic Data Commons using `TCGAbiolinks`.

Clinical TNBC was defined as ER− / PR− / HER2−.

The resulting dataset contained 1,881 miRNA features across 1,207 samples.

### GEO GSE40049

GSE40049 was used for independent validation of TNBC-associated miRNA expression.

- 24 TNBC samples
- 14 adjacent-normal samples

The sequencing count data were analysed using DESeq2.

### GEO GSE86277

GSE86277 was used to investigate TNBC specificity.

The downloaded series matrix used in this analysis contained:

- 63 TNBC samples
- 9 non-TNBC samples

Because this dataset is microarray-based, differential expression was evaluated using `limma`.

## Bioinformatics Workflow

### 1. TCGA Data Acquisition

TCGA-BRCA miRNA expression data were retrieved using `TCGAbiolinks` and GDC data resources.

### 2. Clinical Classification

Samples were classified using ER, PR, and HER2 clinical information into TNBC, HR+/HER2−, HER2+, and normal groups.

### 3. Differential Expression Analysis

Differential expression was performed using DESeq2.

TNBC was compared with:

1. Normal breast tissue
2. HR+/HER2− breast cancer
3. HER2+ breast cancer

Candidate filtering criteria:

- Adjusted p-value < 0.05
- Absolute log2 fold-change ≥ 1
- Consistent direction of expression

### 4. Candidate Selection

The three TCGA comparisons were integrated to identify miRNAs showing consistent differential expression.

Results:

- 218 significant miRNAs in TNBC vs Normal
- 144 significant miRNAs with |log2FC| ≥ 1 in TNBC vs HR+/HER2−
- 114 significant miRNAs with |log2FC| ≥ 1 in TNBC vs HER2+
- 89 candidates were common to the two subtype comparisons
- 34 candidates showed consistent direction across all three comparisons

### 5. Independent Validation

Candidate miRNAs were evaluated in two independent GEO datasets.

#### GSE40049

Four candidates satisfied the predefined replication criteria:

- miR-1304
- miR-19b-2
- miR-18a
- miR-934

#### GSE86277

The four replicated candidates were evaluated for TNBC-specific expression. Among them, miR-934 demonstrated the strongest evidence of TNBC enrichment.

## Key Findings

### TCGA Discovery

A set of 34 directionally consistent candidate miRNAs was identified across TNBC vs Normal, TNBC vs HR+/HER2−, and TNBC vs HER2+ comparisons.

### Independent Validation

Four candidates were independently replicated in GSE40049:

| Candidate | TCGA TNBC vs Normal | GSE40049 | GSE40049 adjusted p |
|---|---:|---:|---:|
| miR-1304 | +1.54 | ~+2.88 | 4.16 × 10⁻⁶ |
| miR-19b-2 | +1.21 | ~+1.59 | 0.00473 |
| miR-18a | +2.80 | ~+1.35 | 0.0395 |
| **miR-934** | **+3.56** | **~+1.79** | **0.0338** |

### TNBC Specificity

In GSE86277, miR-934 showed significantly higher expression in TNBC compared with non-TNBC samples:

- logFC = +2.873
- Adjusted p = 0.0104

### ROC Analysis

miR-934 demonstrated discriminatory potential for TNBC classification:

- AUC = 0.824
- 95% CI = 0.689–0.959

### Survival Analysis

No significant association between miR-934 expression and overall survival was observed:

- Kaplan–Meier log-rank p = 0.822
- Cox HR = 0.912
- 95% CI = 0.401–2.075
- Cox p = 0.827

Therefore, this analysis does not support miR-934 as a prognostic biomarker.

## Functional Analysis

Experimentally supported miRNA-target interactions were investigated using `multiMiR`.

Shared targets among the prioritized miRNA candidates were analysed using Gene Ontology enrichment.

Significantly enriched biological processes included:

- Negative regulation of gene expression
- Regulatory ncRNA-mediated gene silencing

A separate target-gene network was constructed for the prioritized candidate miR-934.

## Main Candidate: miR-934

The integrated analysis prioritized **hsa-miR-934** based on:

- Differential expression in TCGA
- Consistent direction across TCGA comparisons
- Independent replication in GSE40049
- Significant TNBC enrichment in GSE86277
- ROC AUC of approximately 0.824

The findings support miR-934 as a **candidate TNBC classification/diagnostic biomarker**.

However, miR-934 should not be considered a clinically validated biomarker based on this computational analysis alone.

## Figures

- **Figure 1:** Overall study methodology and analysis workflow.
- **Figure 2:** TCGA differential expression analysis using volcano plots.
- **Figure 3:** Heatmap of top differentially expressed miRNAs in TNBC.
- **Figure 4:** Independent validation of replicated candidate miRNAs in GSE40049.
- **Figure 5:** miR-934 expression in TNBC and non-TNBC samples from GSE86277.
- **Figure 6:** ROC curve showing the discriminatory performance of miR-934.
- **Figure 7:** Kaplan–Meier overall-survival analysis of miR-934 expression.
- **Figure 8:** GO enrichment analysis of shared validated target genes.
- **Figure 9:** Experimentally supported miR-934 target-gene network.

## Software and R Packages

The analysis was performed in R.

Major packages used include:

- `TCGAbiolinks`
- `DESeq2`
- `GEOquery`
- `limma`
- `pheatmap`
- `EnhancedVolcano`
- `pROC`
- `survival`
- `survminer`
- `miRBaseConverter`
- `multiMiR`
- `clusterProfiler`
- `org.Hs.eg.db`
- `igraph`

## Limitations

1. The study is based on retrospective publicly available datasets.
2. Independent validation datasets had relatively small sample sizes.
3. GSE86277 contained only 9 non-TNBC samples in the downloaded series matrix.
4. TCGA sequencing and GEO microarray data were generated using different platforms.
5. Only 24 deaths were available for the TCGA TNBC survival analysis.
6. No experimental validation of miR-934 or its targets was performed.
7. miR-934 has previously been reported in association with TNBC; therefore, this project represents integrative computational identification and independent validation rather than discovery of a novel biomarker.

## Conclusion

This study demonstrates an integrative workflow for identifying and independently evaluating miRNA biomarkers associated with TNBC using public transcriptomic datasets.

The analysis prioritized **miR-934**, which showed reproducible TNBC-associated expression across TCGA and GEO datasets and demonstrated discriminatory potential in an independent cohort.

Although miR-934 was not significantly associated with overall survival, its expression pattern supports further investigation as a candidate biomarker for **TNBC classification/diagnostic discrimination**.

Further validation using larger independent cohorts and experimental studies will be required to establish its clinical utility.

## Repository Structure

```text
tnbc-mirna-biomarker-analysis/
├── README.md
├── R/
│   ├── 01_TCGA_Discovery.R
│   ├── 02_GEO_Validation.R
│   ├── 03_Survival_Analysis.R
│   └── 04_miR934_Target_Analysis.R
├── Figures/
├── Results/
└── Docs/
```

## Author

**Aarkin Nanda**  
M.Sc. Biotechnology

## Disclaimer

This project is intended for research and educational purposes. The computational findings do not constitute clinical diagnostic validation or a clinically approved biomarker.
