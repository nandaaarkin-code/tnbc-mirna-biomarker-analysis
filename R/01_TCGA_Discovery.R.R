# ============================================================
# TNBC miRNA Biomarker Discovery - TCGA-BRCA
# Complete analysis script through functional enrichment
#
# Project aim:
# Identify miRNAs differentially expressed in clinically defined
# triple-negative breast cancer (TNBC), compare TNBC with normal
# breast tissue and other breast-cancer receptor groups, identify
# consistent candidate biomarkers, map precursor miRNAs to mature
# miRNAs, retrieve validated targets, and perform GO enrichment.
#
# Data source: TCGA-BRCA miRNA expression + clinical information
# Main statistical method: DESeq2
# ============================================================

# ------------------------------------------------------------
# 0. LOAD PACKAGES
# ------------------------------------------------------------

library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(miRBaseConverter)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# 1. QUERY TCGA-BRCA miRNA DATA
# ------------------------------------------------------------

query_mirna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling"
)

# ------------------------------------------------------------
# 2. DOWNLOAD DATA
# ------------------------------------------------------------
# API/chunked download was used because the standard archive
# download produced an archive-format error.

GDCdownload(
  query_mirna,
  method = "api",
  files.per.chunk = 20
)

# ------------------------------------------------------------
# 3. PREPARE DATA
# ------------------------------------------------------------

mirna_data <- GDCprepare(query_mirna)

# Check the first column / object structure
names(mirna_data)[1]

# ------------------------------------------------------------
# 4. CREATE RAW miRNA COUNT MATRIX
# ------------------------------------------------------------

count_cols <- grep(
  "^read_count_",
  names(mirna_data),
  value = TRUE
)

mirna_counts <- as.matrix(
  mirna_data[, count_cols]
)

rownames(mirna_counts) <- mirna_data$miRNA_ID

colnames(mirna_counts) <- sub(
  "^read_count_",
  "",
  count_cols
)

# Check dimensions
print(dim(mirna_counts))

# ------------------------------------------------------------
# 5. GET TCGA-BRCA PAM50 SUBTYPE INFORMATION
# ------------------------------------------------------------

brca_subtype <- TCGAquery_subtype(
  tumor = "BRCA"
)

table(
  brca_subtype$BRCA_Subtype_PAM50,
  useNA = "ifany"
)

# ------------------------------------------------------------
# 6. MATCH miRNA SAMPLES TO TCGA PATIENTS
# ------------------------------------------------------------

sample_ids <- colnames(mirna_counts)

patient_ids <- substr(
  sample_ids,
  1,
  12
)

subtype_match <- match(
  patient_ids,
  brca_subtype$patient
)

sample_subtype <- brca_subtype$BRCA_Subtype_PAM50[
  subtype_match
]

table(
  sample_subtype,
  useNA = "ifany"
)

# ------------------------------------------------------------
# 7. DETERMINE TCGA SAMPLE TYPE
# ------------------------------------------------------------
# TCGA barcode positions 14-15:
# 01 = Primary Tumor
# 11 = Solid Tissue Normal

sample_type_code <- substr(
  colnames(mirna_counts),
  14,
  15
)

table(sample_type_code)

# ------------------------------------------------------------
# 8. CREATE INITIAL SAMPLE METADATA
# ------------------------------------------------------------

metadata <- data.frame(
  sample = colnames(mirna_counts),
  patient = substr(colnames(mirna_counts), 1, 12),
  sample_type_code = sample_type_code,
  PAM50 = sample_subtype,
  stringsAsFactors = FALSE
)

table(
  metadata$sample_type_code,
  metadata$PAM50,
  useNA = "ifany"
)

# ------------------------------------------------------------
# 9. DOWNLOAD TCGA-BRCA CLINICAL RECEPTOR INFORMATION
# ------------------------------------------------------------

xena_url <- "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"

download.file(
  xena_url,
  destfile = "BRCA_clinicalMatrix.txt",
  mode = "wb"
)

brca_clinical <- read.delim(
  "BRCA_clinicalMatrix.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 10. EXTRACT ER, PR AND HER2 INFORMATION
# ------------------------------------------------------------

clinical_receptor <- brca_clinical[, c(
  "bcr_patient_barcode",
  "ER_Status_nature2012",
  "PR_Status_nature2012",
  "HER2_Final_Status_nature2012"
)]

clinical_receptor <- clinical_receptor[
  !duplicated(clinical_receptor$bcr_patient_barcode),
]

# ------------------------------------------------------------
# 11. MATCH miRNA SAMPLES TO CLINICAL PATIENTS
# ------------------------------------------------------------

miRNA_patient <- substr(
  colnames(mirna_counts),
  1,
  12
)

clinical_match <- match(
  miRNA_patient,
  clinical_receptor$bcr_patient_barcode
)

miRNA_clinical <- clinical_receptor[
  clinical_match,
]

rownames(miRNA_clinical) <- colnames(mirna_counts)

miRNA_clinical$sample <- colnames(mirna_counts)
miRNA_clinical$patient <- miRNA_patient

# Add PAM50 information
miRNA_clinical$PAM50 <- sample_subtype

# Add TCGA sample type
miRNA_clinical$sample_type_code <- sample_type_code

# ------------------------------------------------------------
# 12. DEFINE CLINICAL BREAST-CANCER GROUPS
# ------------------------------------------------------------
# TNBC = ER-negative + PR-negative + HER2-negative
# HER2+ = HER2-positive
# HR+/HER2- = ER-positive + HER2-negative

miRNA_clinical$clinical_group <- NA_character_

# TNBC
miRNA_clinical$clinical_group[
  miRNA_clinical$ER_Status_nature2012 == "Negative" &
  miRNA_clinical$PR_Status_nature2012 == "Negative" &
  miRNA_clinical$HER2_Final_Status_nature2012 == "Negative"
] <- "TNBC"

# HER2-positive
miRNA_clinical$clinical_group[
  miRNA_clinical$HER2_Final_Status_nature2012 == "Positive"
] <- "HER2+"

# HR-positive / HER2-negative
miRNA_clinical$clinical_group[
  miRNA_clinical$ER_Status_nature2012 == "Positive" &
  miRNA_clinical$HER2_Final_Status_nature2012 == "Negative"
] <- "HR+/HER2-"

# Everything else
miRNA_clinical$clinical_group[
  is.na(miRNA_clinical$clinical_group)
] <- "Other/Unknown"

# Check group sizes
table(miRNA_clinical$clinical_group)

# Compare clinical TNBC with PAM50
# Clinical TNBC and PAM50 Basal are not identical.
table(
  miRNA_clinical$clinical_group,
  miRNA_clinical$PAM50,
  useNA = "ifany"
)

# ------------------------------------------------------------
# 13. PREPARE PRIMARY-TUMOR CLINICAL GROUPS FOR DESeq2
# ------------------------------------------------------------

# Keep primary tumor samples from the three clinically defined groups.
tumor_clinical_metadata <- miRNA_clinical[
  miRNA_clinical$sample_type_code == "01" &
  miRNA_clinical$clinical_group %in%
    c("TNBC", "HR+/HER2-", "HER2+"),
]

# Keep one primary tumor sample per patient.
tumor_clinical_metadata <- tumor_clinical_metadata[
  !duplicated(tumor_clinical_metadata$patient),
]

table(tumor_clinical_metadata$clinical_group)

# ------------------------------------------------------------
# 14. ALIGN COUNT MATRIX WITH CLINICAL METADATA
# ------------------------------------------------------------

clinical_samples <- tumor_clinical_metadata$sample

mirna_counts_clinical <- mirna_counts[
  ,
  clinical_samples,
  drop = FALSE
]

# Confirm exact alignment
print(
  all(
    colnames(mirna_counts_clinical) ==
      tumor_clinical_metadata$sample
  )
)

# ------------------------------------------------------------
# 15. CREATE DESeq2 DATASET FOR THREE CLINICAL GROUPS
# ------------------------------------------------------------

clinical_coldata <- tumor_clinical_metadata

rownames(clinical_coldata) <-
  clinical_coldata$sample

clinical_coldata$clinical_group <- factor(
  clinical_coldata$clinical_group,
  levels = c(
    "HR+/HER2-",
    "HER2+",
    "TNBC"
  )
)

dds_clinical <- DESeqDataSetFromMatrix(
  countData = round(mirna_counts_clinical),
  colData = clinical_coldata,
  design = ~ clinical_group
)

# ------------------------------------------------------------
# 16. FILTER LOW-COUNT miRNAs
# ------------------------------------------------------------

# Keep miRNAs with at least 10 total reads across the dataset.
dds_clinical <- dds_clinical[
  rowSums(counts(dds_clinical)) >= 10,
]

print(nrow(dds_clinical))

# ------------------------------------------------------------
# 17. RUN DESeq2
# ------------------------------------------------------------

dds_clinical <- DESeq(dds_clinical)

# ============================================================
# 18. DIFFERENTIAL EXPRESSION: TNBC vs HR+/HER2-
# ============================================================

res_tnbc_hr <- results(
  dds_clinical,
  contrast = c(
    "clinical_group",
    "TNBC",
    "HR+/HER2-"
  )
)

res_tnbc_hr_df <- as.data.frame(
  res_tnbc_hr
)

res_tnbc_hr_df$miRNA <-
  rownames(res_tnbc_hr_df)

# Significant miRNAs
sig_tnbc_hr <- res_tnbc_hr_df[
  !is.na(res_tnbc_hr_df$padj) &
  res_tnbc_hr_df$padj < 0.05 &
  abs(res_tnbc_hr_df$log2FoldChange) >= 1,
]

print(nrow(sig_tnbc_hr))

# ============================================================
# 19. DIFFERENTIAL EXPRESSION: TNBC vs HER2+
# ============================================================

res_tnbc_her2 <- results(
  dds_clinical,
  contrast = c(
    "clinical_group",
    "TNBC",
    "HER2+"
  )
)

res_tnbc_her2_df <- as.data.frame(
  res_tnbc_her2
)

res_tnbc_her2_df$miRNA <-
  rownames(res_tnbc_her2_df)

# Significant miRNAs
sig_tnbc_her2 <- res_tnbc_her2_df[
  !is.na(res_tnbc_her2_df$padj) &
  res_tnbc_her2_df$padj < 0.05 &
  abs(res_tnbc_her2_df$log2FoldChange) >= 1,
]

print(nrow(sig_tnbc_her2))

# ============================================================
# 20. PREPARE TNBC AND NORMAL SAMPLES
# ============================================================

# Select primary TNBC tumors.
tnbc_metadata <- miRNA_clinical[
  miRNA_clinical$sample_type_code == "01" &
  miRNA_clinical$clinical_group == "TNBC",
]

# One sample per patient.
tnbc_metadata <- tnbc_metadata[
  !duplicated(tnbc_metadata$patient),
]

# Select solid tissue normal samples.
normal_metadata <- miRNA_clinical[
  miRNA_clinical$sample_type_code == "11",
]

# One normal sample per patient.
normal_metadata <- normal_metadata[
  !duplicated(normal_metadata$patient),
]

# ------------------------------------------------------------
# 21. BUILD TNBC VS NORMAL METADATA
# ------------------------------------------------------------

tnbc_metadata$condition <- "TNBC"
normal_metadata$condition <- "Normal"

tnbc_normal_coldata <- rbind(
  tnbc_metadata,
  normal_metadata
)

tnbc_normal_coldata$condition <- factor(
  tnbc_normal_coldata$condition,
  levels = c("Normal", "TNBC")
)

rownames(tnbc_normal_coldata) <-
  tnbc_normal_coldata$sample

table(tnbc_normal_coldata$condition)

# ------------------------------------------------------------
# 22. ALIGN TNBC VS NORMAL COUNT MATRIX
# ------------------------------------------------------------

tnbc_normal_samples <-
  tnbc_normal_coldata$sample

mirna_counts_tnbc_normal <- mirna_counts[
  ,
  tnbc_normal_samples,
  drop = FALSE
]

print(
  all(
    colnames(mirna_counts_tnbc_normal) ==
      tnbc_normal_coldata$sample
  )
)

# ------------------------------------------------------------
# 23. CREATE TNBC VS NORMAL DESeq2 DATASET
# ------------------------------------------------------------

dds_tnbc_normal <- DESeqDataSetFromMatrix(
  countData = round(mirna_counts_tnbc_normal),
  colData = tnbc_normal_coldata,
  design = ~ condition
)

# Filter low-count miRNAs
dds_tnbc_normal <- dds_tnbc_normal[
  rowSums(counts(dds_tnbc_normal)) >= 10,
]

print(nrow(dds_tnbc_normal))

# ------------------------------------------------------------
# 24. RUN TNBC VS NORMAL DESeq2
# ------------------------------------------------------------

dds_tnbc_normal <- DESeq(
  dds_tnbc_normal
)

res_tnbc_normal <- results(
  dds_tnbc_normal,
  contrast = c(
    "condition",
    "TNBC",
    "Normal"
  )
)

res_tnbc_normal_df <- as.data.frame(
  res_tnbc_normal
)

res_tnbc_normal_df$miRNA <-
  rownames(res_tnbc_normal_df)

# Significant miRNAs
sig_tnbc_normal <- res_tnbc_normal_df[
  !is.na(res_tnbc_normal_df$padj) &
  res_tnbc_normal_df$padj < 0.05 &
  abs(res_tnbc_normal_df$log2FoldChange) >= 1,
]

print(nrow(sig_tnbc_normal))

# ============================================================
# 25. FIND miRNAs SIGNIFICANT IN ALL THREE COMPARISONS
# ============================================================

common_sig <- Reduce(
  intersect,
  list(
    sig_tnbc_normal$miRNA,
    sig_tnbc_hr$miRNA,
    sig_tnbc_her2$miRNA
  )
)

print(length(common_sig))

# Create a combined table
common_table <- data.frame(
  miRNA = common_sig,
  stringsAsFactors = FALSE
)

common_table$log2FC_TNBC_vs_Normal <-
  res_tnbc_normal_df$log2FoldChange[
    match(
      common_table$miRNA,
      res_tnbc_normal_df$miRNA
    )
  ]

common_table$log2FC_TNBC_vs_HR <-
  res_tnbc_hr_df$log2FoldChange[
    match(
      common_table$miRNA,
      res_tnbc_hr_df$miRNA
    )
  ]

common_table$log2FC_TNBC_vs_HER2 <-
  res_tnbc_her2_df$log2FoldChange[
    match(
      common_table$miRNA,
      res_tnbc_her2_df$miRNA
    )
  ]

# ------------------------------------------------------------
# 26. KEEP ONLY DIRECTION-CONSISTENT CANDIDATES
# ------------------------------------------------------------

common_table$direction_consistent <-
  (
    sign(common_table$log2FC_TNBC_vs_Normal) ==
    sign(common_table$log2FC_TNBC_vs_HR)
  ) &
  (
    sign(common_table$log2FC_TNBC_vs_Normal) ==
    sign(common_table$log2FC_TNBC_vs_HER2)
  )

common_consistent <- common_table[
  common_table$direction_consistent,
]

print(nrow(common_consistent))

# ------------------------------------------------------------
# 27. RANK CONSISTENT CANDIDATES BY EFFECT SIZE
# ------------------------------------------------------------

common_consistent$min_abs_log2FC <- apply(
  abs(
    common_consistent[
      ,
      c(
        "log2FC_TNBC_vs_Normal",
        "log2FC_TNBC_vs_HR",
        "log2FC_TNBC_vs_HER2"
      )
    ]
  ),
  1,
  min
)

common_consistent$mean_abs_log2FC <- rowMeans(
  abs(
    common_consistent[
      ,
      c(
        "log2FC_TNBC_vs_Normal",
        "log2FC_TNBC_vs_HR",
        "log2FC_TNBC_vs_HER2"
      )
    ]
  )
)

common_consistent <- common_consistent[
  order(
    -common_consistent$min_abs_log2FC
  ),
]

print(common_consistent)

# ============================================================
# 28. SAVE DESeq2 RESULT TABLES
# ============================================================

write.csv(
  res_tnbc_normal_df,
  "DESeq2_TNBC_vs_Normal_miRNA.csv",
  row.names = FALSE
)

write.csv(
  res_tnbc_hr_df,
  "DESeq2_TNBC_vs_HR_HER2neg_miRNA.csv",
  row.names = FALSE
)

write.csv(
  res_tnbc_her2_df,
  "DESeq2_TNBC_vs_HER2pos_miRNA.csv",
  row.names = FALSE
)

write.csv(
  common_consistent,
  "TNBC_consistent_miRNA_candidates.csv",
  row.names = FALSE
)

# ============================================================
# 29. VOLCANO PLOT: TNBC VS NORMAL
# ============================================================

EnhancedVolcano(
  res_tnbc_normal_df,
  lab = res_tnbc_normal_df$miRNA,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "TNBC vs Normal",
  subtitle = "TCGA-BRCA miRNA",
  legendPosition = "right"
)

# ============================================================
# 30. VOLCANO PLOT: TNBC VS HR+/HER2-
# ============================================================

EnhancedVolcano(
  res_tnbc_hr_df,
  lab = res_tnbc_hr_df$miRNA,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "TNBC vs HR+/HER2-",
  subtitle = "TCGA-BRCA miRNA",
  legendPosition = "right"
)

# ============================================================
# 31. VOLCANO PLOT: TNBC VS HER2+
# ============================================================

EnhancedVolcano(
  res_tnbc_her2_df,
  lab = res_tnbc_her2_df$miRNA,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "TNBC vs HER2+",
  subtitle = "TCGA-BRCA miRNA",
  legendPosition = "right"
)

# ============================================================
# 32. HEATMAP: TOP 30 TNBC VS NORMAL miRNAs
# ============================================================

# The default vst() failed because the filtered dataset contained
# fewer rows than the default nsub. Use the full transformation instead.
vsd_tnbc_normal <- varianceStabilizingTransformation(
  dds_tnbc_normal,
  blind = FALSE
)

top30_mirnas <- head(
  sig_tnbc_normal[
    order(sig_tnbc_normal$padj),
  ],
  30
)

heatmap_mat <- assay(
  vsd_tnbc_normal
)[
  top30_mirnas$miRNA,
  ,
  drop = FALSE
]

# Scale each miRNA by row
heatmap_mat_scaled <- t(
  scale(t(heatmap_mat))
)

pheatmap(
  heatmap_mat_scaled,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 9,
  main = "Top 30 Differentially Expressed miRNAs\nTNBC vs Normal"
)

# ============================================================
# 33. MAP TCGA PRECURSOR miRNA IDs TO MATURE miRNAs
# ============================================================

candidate_ids <- common_consistent$miRNA

candidate_mature <- miRNA_PrecursorToMature(
  candidate_ids
)

print(candidate_mature)

# ============================================================
# 34. CREATE INTEGRATED CANDIDATE TABLE
# ============================================================

candidate_table <- candidate_mature

candidate_table$log2FC_TNBC_vs_Normal <-
  common_consistent$log2FC_TNBC_vs_Normal[
    match(
      candidate_table$OriginalName,
      common_consistent$miRNA
    )
  ]

candidate_table$log2FC_TNBC_vs_HR <-
  common_consistent$log2FC_TNBC_vs_HR[
    match(
      candidate_table$OriginalName,
      common_consistent$miRNA
    )
  ]

candidate_table$log2FC_TNBC_vs_HER2 <-
  common_consistent$log2FC_TNBC_vs_HER2[
    match(
      candidate_table$OriginalName,
      common_consistent$miRNA
    )
  ]

candidate_table$min_abs_log2FC <-
  common_consistent$min_abs_log2FC[
    match(
      candidate_table$OriginalName,
      common_consistent$miRNA
    )
  ]

candidate_table$mean_abs_log2FC <-
  common_consistent$mean_abs_log2FC[
    match(
      candidate_table$OriginalName,
      common_consistent$miRNA
    )
  ]

write.csv(
  candidate_table,
  "TNBC_miRNA_candidate_mature_mapping.csv",
  row.names = FALSE
)

# ============================================================
# 35. SELECT TOP MATURE miRNA CANDIDATES FOR TARGET ANALYSIS
# ============================================================

# These are the strongest preliminary candidates based on the
# integrated TCGA differential-expression results.
# Mature arm selection will be refined during later validation.

top_precursors <- c(
  "hsa-mir-577",
  "hsa-mir-522",
  "hsa-mir-934",
  "hsa-mir-516a-1",
  "hsa-mir-516a-2",
  "hsa-mir-138-1",
  "hsa-mir-138-2"
)

top_mature_map <- candidate_mature[
  candidate_mature$OriginalName %in% top_precursors,
]

print(top_mature_map)

# Extract unique mature IDs from the mapping.
top_mature_ids <- unique(
  c(
    top_mature_map$Mature1,
    top_mature_map$Mature2
  )
)

top_mature_ids <- top_mature_ids[
  !is.na(top_mature_ids) &
  top_mature_ids != ""
]

print(top_mature_ids)

# ============================================================
# 36. RETRIEVE EXPERIMENTALLY VALIDATED miRNA TARGETS
# ============================================================

validated_top <- get_multimir(
  org = "hsa",
  mirna = top_mature_ids,
  table = "validated",
  summary = TRUE
)

# Extract validated interaction table
validated_df <- validated_top@data

validated_df <- validated_df[
  validated_df$type == "validated",
]

# ============================================================
# 37. UNIQUE VALIDATED TARGET GENES
# ============================================================

target_genes <- unique(
  validated_df$target_symbol
)

print(length(target_genes))

# ============================================================
# 38. SHARED TARGET GENES
# ============================================================
# Genes targeted by at least two candidate miRNAs.

target_counts <- table(
  validated_df$target_symbol
)

shared_targets <- names(
  target_counts[
    target_counts >= 2
  ]
)

print(length(shared_targets))

# ============================================================
# 39. TARGET COUNT SUMMARY PER miRNA
# ============================================================

target_summary <- aggregate(
  target_symbol ~ mature_mirna_id,
  data = validated_df,
  FUN = function(x) length(unique(x))
)

target_summary <- target_summary[
  order(-target_summary$target_symbol),
]

print(target_summary)

write.csv(
  target_summary,
  "TNBC_candidate_miRNA_validated_target_summary.csv",
  row.names = FALSE
)

# ============================================================
# 40. CONVERT SHARED TARGET SYMBOLS TO ENTREZ IDs
# ============================================================

target_conversion <- bitr(
  shared_targets,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

target_conversion <- target_conversion[
  !duplicated(target_conversion$ENTREZID),
]

print(nrow(target_conversion))

# ============================================================
# 41. CREATE APPROPRIATE GO BACKGROUND
# ============================================================
# Use all experimentally validated targets of the selected
# candidate miRNAs as the background rather than the whole genome.

background_genes <- unique(
  validated_df$target_symbol
)

background_conversion <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_conversion <- background_conversion[
  !duplicated(background_conversion$ENTREZID),
]

print(nrow(background_conversion))

# ============================================================
# 42. GO BIOLOGICAL PROCESS ENRICHMENT
# ============================================================

ego_bp_bg <- enrichGO(
  gene = target_conversion$ENTREZID,
  universe = background_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

print(ego_bp_bg)

# View enrichment results
as.data.frame(ego_bp_bg)

# ============================================================
# 43. GO DOT PLOT
# ============================================================

dotplot(
  ego_bp_bg,
  showCategory = 10,
  title = "GO Biological Process Enrichment of Shared miRNA Targets"
)

# ============================================================
# 44. SAVE GO RESULTS
# ============================================================

write.csv(
  as.data.frame(ego_bp_bg),
  "GO_BP_shared_miRNA_targets.csv",
  row.names = FALSE
)

# ============================================================
# 45. OPTIONAL: SAVE ALL OBJECTS FOR REPRODUCIBILITY
# ============================================================

save(
  mirna_counts,
  miRNA_clinical,
  dds_clinical,
  dds_tnbc_normal,
  res_tnbc_normal_df,
  res_tnbc_hr_df,
  res_tnbc_her2_df,
  common_consistent,
  candidate_table,
  validated_df,
  ego_bp_bg,
  file = "TNBC_miRNA_TCGA_analysis_objects.RData"
)

# ============================================================
# END OF CURRENT TCGA ANALYSIS
# ============================================================
# Next planned stage:
# Independent GEO validation of the strongest candidate miRNAs.
# Do NOT label candidates as validated TNBC biomarkers until
# independent validation has been completed.
# ============================================================
