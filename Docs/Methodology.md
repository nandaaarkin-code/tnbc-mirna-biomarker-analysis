# Methodology

## Study Design

An integrative bioinformatics workflow was used to identify and independently validate candidate miRNA biomarkers associated with triple-negative breast cancer (TNBC).

The analysis consisted of four major stages:

**TCGA discovery → candidate selection → GEO validation → functional and clinical evaluation**

## 1. TCGA-BRCA miRNA Discovery

miRNA expression data for TCGA-BRCA were obtained from the Genomic Data Commons using the `TCGAbiolinks` package in R.

The dataset contained 1,881 miRNA features across 1,207 samples.

Clinical information was obtained from the TCGA-BRCA clinical matrix. Tumor samples were classified using clinical ER, PR and HER2 status:

- **TNBC:** ER-negative, PR-negative and HER2-negative
- **HR+/HER2−:** ER-positive and HER2-negative
- **HER2+:** HER2-positive

Normal breast samples were also retained for comparison.

Importantly, clinical TNBC classification was used rather than relying solely on the PAM50 Basal subtype.

## 2. Differential miRNA Expression Analysis

Differential expression analysis was performed using **DESeq2**.

Three comparisons were performed:

1. TNBC vs Normal
2. TNBC vs HR+/HER2−
3. TNBC vs HER2+

MiRNAs were considered significant when they satisfied:

- Adjusted p-value < 0.05
- Absolute log2 fold change ≥ 1

Candidate miRNAs were further prioritized based on:

- Significance across multiple comparisons
- Sufficient effect size
- Consistent direction of expression change

Precursor miRNA identifiers were converted to mature miRNA identifiers using `miRBaseConverter` where applicable.

## 3. Independent GEO Validation

Two independent GEO datasets were used for validation.

### GSE40049

GSE40049 contains small RNA sequencing data from TNBC and adjacent normal breast tissue.

The processed count files were downloaded from GEO and combined into a miRNA count matrix. Differential expression was performed using DESeq2.

Candidates were considered replicated when they showed:

- Adjusted p-value < 0.05
- Absolute log2 fold change ≥ 1
- The same direction of expression change as observed in TCGA

### GSE86277

GSE86277 contains miRNA expression data generated using a microarray platform.

Because these data are normalized microarray expression values rather than sequencing counts, differential expression was evaluated using the **limma** framework.

TNBC samples were compared with the available non-TNBC breast cancer subtypes.

This analysis was used to assess whether the replicated candidates showed TNBC-associated expression patterns independent of the TCGA dataset and sequencing platform.

## 4. Candidate Prioritization

Candidates that demonstrated reproducibility between TCGA and GSE40049 were examined further in GSE86277.

Among the replicated candidates, **hsa-miR-934** showed significant enrichment in TNBC compared with non-TNBC samples in GSE86277.

Therefore, hsa-miR-934 was selected as the lead candidate for downstream evaluation.

## 5. ROC Analysis

Receiver operating characteristic (ROC) analysis was performed using the `pROC` package to evaluate the ability of miR-934 expression to distinguish TNBC from non-TNBC samples in GSE86277.

The area under the ROC curve (AUC), confidence interval, sensitivity and specificity were calculated.

## 6. Survival Analysis

Overall survival analysis was performed in the TCGA TNBC cohort.

Patients were divided into high- and low-miR-934 expression groups using the median expression as the cutoff.

Kaplan-Meier survival analysis and Cox proportional hazards regression were performed.

The proportional hazards assumption was also evaluated.

## 7. Target and Functional Analysis

Validated miR-934 target interactions were retrieved using the `multiMiR` package.

Validated target genes were subjected to functional enrichment analysis using `clusterProfiler`.

Gene Ontology biological-process enrichment was examined to investigate biological processes associated with the validated target set.

A target interaction network was also constructed using `igraph` to visualize selected validated miR-934 target genes.

## 8. Software and R Packages

The analysis was performed in R using packages including:

- `TCGAbiolinks`
- `DESeq2`
- `GEOquery`
- `limma`
- `miRBaseConverter`
- `multiMiR`
- `clusterProfiler`
- `org.Hs.eg.db`
- `pROC`
- `survival`
- `survminer`
- `pheatmap`
- `EnhancedVolcano`
- `igraph`

## Overall Workflow

```text
TCGA-BRCA miRNA data
        ↓
Clinical TNBC classification
        ↓
TNBC vs Normal
TNBC vs HR+/HER2−
TNBC vs HER2+
        ↓
Differential miRNA analysis
        ↓
Robust candidate selection
        ↓
GSE40049 independent validation
        ↓
GSE86277 TNBC specificity analysis
        ↓
Lead candidate: hsa-miR-934
        ↓
ROC / AUC analysis
        ↓
Overall survival analysis
        ↓
Validated target analysis
        ↓
GO enrichment + target network
