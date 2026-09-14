# Results

## 1. TCGA-BRCA miRNA Dataset

TCGA-BRCA miRNA expression data contained **1,881 miRNA features across 1,207 samples**.

Clinical classification identified:

- TNBC: 133 samples
- HR+/HER2−: 542 samples
- HER2+: 125 samples
- Other/Unknown: 407 samples

For the differential-expression analyses, after restricting to primary tumor samples and removing duplicate patients, the final comparison cohorts contained **119 TNBC, 473 HR+/HER2− and 111 HER2+ samples**.

Normal breast samples were analyzed separately.

## 2. Differential Expression in TCGA

Three subtype comparisons were performed using DESeq2:

- TNBC vs Normal
- TNBC vs HR+/HER2−
- TNBC vs HER2+

Using adjusted p-value < 0.05 and |log2FC| ≥ 1:

| Comparison | Significant miRNAs |
|---|---:|
| TNBC vs Normal | 218 |
| TNBC vs HR+/HER2− | 144 |
| TNBC vs HER2+ | 114 |

A total of **89 miRNAs** were significantly differentially expressed in both TNBC vs HR+/HER2− and TNBC vs HER2+ comparisons.

Of these, **88 showed consistent direction of expression**, indicating that the majority of these candidates demonstrated a reproducible TNBC-associated expression pattern across the two breast cancer subtype comparisons.

When TNBC vs Normal was also incorporated, **38 miRNAs** were significant across all three comparisons.

Among these, **34 showed consistent direction of expression across all three comparisons** and were retained as robust TCGA candidate miRNAs.

## 3. Prioritized TCGA Candidates

Candidates were ranked using their minimum absolute log2 fold change across the three TCGA comparisons.

Several candidates showed strong and consistent expression changes, including:

- hsa-mir-577
- hsa-mir-522
- hsa-mir-516a
- hsa-mir-934
- hsa-mir-519a
- hsa-mir-138
- hsa-mir-137
- hsa-mir-551a
- hsa-mir-135b

The precursor identifiers were subsequently mapped to mature miRNA identifiers where appropriate.

## 4. Independent Validation Using GSE40049

GSE40049 was used as an independent TNBC-versus-normal small RNA sequencing dataset.

The dataset contained:

- 24 TNBC samples
- 14 adjacent normal samples

DESeq2 identified **209 significantly differentially expressed miRNAs** at adjusted p-value < 0.05.

Candidate replication required:

- adjusted p-value < 0.05
- |log2FC| ≥ 1
- same direction of expression as the TCGA TNBC-versus-normal comparison

Four TCGA candidates satisfied these criteria:

| miRNA | TCGA log2FC | GEO log2FC | GEO adjusted p-value |
|---|---:|---:|---:|
| hsa-mir-1304 | +1.54 | ~+2.88 | 4.16 × 10⁻⁶ |
| hsa-mir-19b-2 | +1.21 | ~+1.59 | 0.00473 |
| hsa-mir-18a | +2.80 | ~+1.35 | 0.0395 |
| hsa-mir-934 | +3.56 | ~+1.79 | 0.0338 |

These four miRNAs therefore showed independent replication of the TNBC-associated expression pattern.

## 5. TNBC Specificity in GSE86277

GSE86277 was used to investigate whether the replicated candidates were preferentially associated with TNBC relative to other breast cancer subtypes.

The downloaded series matrix contained:

- 63 TNBC samples
- 9 non-TNBC samples

Because GSE86277 contains normalized microarray expression data, the analysis was performed using the limma framework.

Among the four replicated candidates, **hsa-miR-934** showed the strongest evidence of TNBC-associated expression:

| Candidate | logFC | Adjusted p-value |
|---|---:|---:|
| hsa-miR-934 | +2.873 | 0.0104 |
| hsa-miR-18a | +1.10 | 0.0867 |
| hsa-miR-1304 | +0.24 | 0.2600 |
| hsa-mir-19b-2 | ~0 | ~0.991 |

Thus, hsa-miR-934 was selected as the lead candidate for downstream evaluation.

## 6. Diagnostic Classification Performance of miR-934

ROC analysis was performed using GSE86277.

The area under the ROC curve was:

**AUC = 0.824**

with a 95% confidence interval of:

**0.689–0.959**

The optimal Youden threshold produced:

- Sensitivity: 52.4%
- Specificity: 100%

The AUC indicates that miR-934 expression showed good discriminatory ability between TNBC and the available non-TNBC samples in this dataset.

However, the non-TNBC group was small, so this result should be considered supportive rather than definitive clinical validation.

## 7. Overall Survival Analysis

Overall survival analysis was performed in the TCGA TNBC cohort.

The final survival dataset contained:

- 135 TNBC patients
- 24 observed deaths
- 111 censored observations

Patients were divided into high- and low-miR-934 expression groups using the median expression cutoff.

Kaplan-Meier analysis showed no statistically significant difference in overall survival:

**Log-rank p = 0.822**

Cox proportional hazards analysis gave:

**HR = 0.912**

95% CI:

**0.401–2.075**

Cox p-value:

**0.827**

The proportional hazards assumption was not violated (p = 0.74).

Therefore, miR-934 was **not supported as a prognostic overall-survival biomarker** in this analysis.

## 8. Validated Target Analysis

Validated miR-934 target interactions were retrieved using multiMiR.

The analysis identified **55 unique validated target genes** for hsa-miR-934 in the queried database.

Functional enrichment analysis was subsequently performed to investigate biological processes associated with the validated targets.

The broader shared-target analysis identified enrichment of:

1. Negative regulation of gene expression
2. Regulatory ncRNA-mediated gene silencing

Both terms remained significant after multiple-testing correction.

KEGG analysis did not identify pathways that remained statistically significant after FDR correction.

A miR-934 target interaction network was also generated to visualize selected validated target genes.

## 9. Overall Candidate Selection

The integrated analysis followed the evidence progression:

**TCGA discovery**

↓  

**38 candidates significant across TNBC vs Normal, HR+/HER2− and HER2+**

↓  

**34 candidates with consistent direction**

↓  

**GSE40049 replication**

↓  

**4 replicated candidates**

↓  

**GSE86277 TNBC specificity**

↓  

**hsa-miR-934 prioritized**

↓  

**ROC analysis: AUC 0.824**

↓  

**Survival analysis: no significant association**

Therefore, hsa-miR-934 was identified as a **candidate TNBC classification/diagnostic biomarker**, rather than a prognostic biomarker.

## 10. Figures

The major findings are represented by the following figures:

- **Figure 1:** Methodology workflow
- **Figure 2:** TCGA differential-expression volcano plots
- **Figure 3:** TCGA miRNA heatmap
- **Figure 4:** GSE40049 candidate validation
- **Figure 5:** miR-934 TNBC specificity analysis
- **Figure 6:** miR-934 ROC curve
- **Figure 7:** miR-934 overall survival analysis
- **Figure 8:** GO enrichment of shared validated targets
- **Figure 9:** miR-934 validated target network

## 11. Overall Interpretation

The integrative analysis identified a set of miRNAs consistently associated with TNBC across TCGA comparisons. Independent validation reduced this candidate set to four miRNAs, of which hsa-miR-934 demonstrated the strongest evidence for TNBC-associated expression in an additional independent dataset.

The ROC analysis further supported its potential for TNBC classification, whereas survival analysis did not support a prognostic role.

Overall, hsa-miR-934 represents a **computationally supported candidate TNBC biomarker requiring further validation in larger and clinically independent cohorts and experimental studies**.
