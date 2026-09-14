# ============================================================
# GEO VALIDATION ANALYSIS
# TNBC miRNA Biomarker Project
#
# Datasets:
#   GSE40049  - TNBC vs adjacent Normal
#   GSE86277  - TNBC vs Non-TNBC
#
# Working GEO code used so far
# ============================================================

# 1. GSE40049
library(GEOquery)

gse40049 <- getGEO("GSE40049", GSEMatrix = TRUE)

length(gse40049)
dim(exprs(gse40049[[1]]))

pData(gse40049[[1]])[
  , c("title", "geo_accession", "source_name_ch1")
]

# 2. Download supplementary processed files
getGEOSuppFiles("GSE40049", makeDirectory = TRUE)

untar(
  "GSE40049/GSE40049_RAW.tar",
  exdir = "GSE40049/raw"
)

length(list.files("GSE40049/raw"))
head(list.files("GSE40049/raw"))

# 3. Inspect one processed file
test_file <- list.files(
  "GSE40049/raw",
  full.names = TRUE
)[1]

test_data <- read.delim(
  gzfile(test_file),
  header = TRUE,
  stringsAsFactors = FALSE
)

head(test_data)
dim(test_data)
colnames(test_data)

# 4. Extract miRNA IDs and counts
test_processed <- data.frame(
  miRNA = sub(
    ".*ID=([^;]+);.*",
    "\\1",
    test_data[[9]]
  ),
  count = as.numeric(test_data[[8]])
)

head(test_processed)
summary(test_processed$count)

# 5. Collapse duplicate miRNA entries
test_counts <- aggregate(
  count ~ miRNA,
  data = test_processed,
  FUN = sum
)

head(test_counts)
nrow(test_counts)

head(
  test_counts[
    order(test_counts$count, decreasing = TRUE),
  ],
  10
)

# 6. Get all GSE40049 files
geo_files <- list.files(
  "GSE40049/raw",
  pattern = "\\.txt\\.gz$",
  full.names = TRUE
)

length(geo_files)

# 7. Function to process one sample
read_geo_miRNA <- function(file) {

  dat <- read.delim(
    gzfile(file),
    header = TRUE,
    stringsAsFactors = FALSE
  )

  result <- data.frame(
    miRNA = sub(
      ".*ID=([^;]+);.*",
      "\\1",
      dat[[9]]
    ),
    count = as.numeric(dat[[8]])
  )

  result <- result[
    !is.na(result$miRNA) &
    !is.na(result$count),
  ]

  result <- aggregate(
    count ~ miRNA,
    data = result,
    FUN = sum
  )

  return(result)
}

# 8. Process all samples
geo_list <- lapply(
  geo_files,
  read_geo_miRNA
)

names(geo_list) <- sub(
  "_.*",
  "",
  basename(geo_files)
)

length(geo_list)

# 9. Create complete count matrix
all_miRNAs <- sort(
  unique(
    unlist(
      lapply(
        geo_list,
        function(x) x$miRNA
      )
    )
  )
)

geo_counts <- matrix(
  0,
  nrow = length(all_miRNAs),
  ncol = length(geo_list),
  dimnames = list(
    all_miRNAs,
    names(geo_list)
  )
)

for (i in seq_along(geo_list)) {

  x <- geo_list[[i]]

  geo_counts[
    x$miRNA,
    i
  ] <- x$count
}

dim(geo_counts)
head(geo_counts[, 1:5])

# 10. TNBC / Normal metadata
geo_metadata <- pData(
  gse40049[[1]]
)

geo_metadata$group <- ifelse(
  grepl("_T$", geo_metadata$title),
  "TNBC",
  "Normal"
)

table(geo_metadata$group)

# 11. Align metadata
all(
  colnames(geo_counts) %in%
    geo_metadata$geo_accession
)

geo_metadata <- geo_metadata[
  match(
    colnames(geo_counts),
    geo_metadata$geo_accession
  ),
]

all(
  colnames(geo_counts) ==
    geo_metadata$geo_accession
)

table(geo_metadata$group)

# 12. Identify TCGA candidates in GSE40049
# candidate_ids must already exist from the TCGA analysis.

id_check <- intersect(
  rownames(geo_counts),
  candidate_ids
)

length(id_check)
id_check

# 13. Extract candidate counts
geo_candidate_counts <- geo_counts[
  id_check,
]

dim(geo_candidate_counts)

summary(
  rowSums(geo_candidate_counts)
)

head(
  geo_candidate_counts[
    order(
      rowSums(geo_candidate_counts),
      decreasing = TRUE
    ),
  ],
  10
)

# 14. Prepare counts for DESeq2
library(DESeq2)

geo_counts <- round(geo_counts)

keep_geo <- rowSums(geo_counts) >= 10

geo_counts_filtered <- geo_counts[
  keep_geo,
]

dim(geo_counts_filtered)

# 15. DESeq2 metadata
geo_coldata <- data.frame(
  row.names = geo_metadata$geo_accession,
  group = factor(
    geo_metadata$group,
    levels = c("Normal", "TNBC")
  )
)

all(
  colnames(geo_counts_filtered) ==
    rownames(geo_coldata)
)

# 16. Run DESeq2
dds_geo <- DESeqDataSetFromMatrix(
  countData = geo_counts_filtered,
  colData = geo_coldata,
  design = ~ group
)

dds_geo

dds_geo <- DESeq(dds_geo)

res_geo <- results(
  dds_geo,
  contrast = c(
    "group",
    "TNBC",
    "Normal"
  )
)

res_geo <- as.data.frame(res_geo)

res_geo <- res_geo[
  order(res_geo$padj),
]

head(res_geo, 10)

sum(
  res_geo$padj < 0.05,
  na.rm = TRUE
)

# 17. Validate TCGA candidates
# res_TNBC_vs_Normal must already exist from TCGA analysis.

tcga_tnbc_normal <- as.data.frame(
  res_TNBC_vs_Normal
)

geo_validation <- res_geo[
  rownames(res_geo) %in% id_check,
  c(
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj"
  )
]

geo_validation$TCGA_log2FC <-
  tcga_tnbc_normal[
    match(
      rownames(geo_validation),
      rownames(tcga_tnbc_normal)
    ),
    "log2FoldChange"
  ]

geo_validation$GEO_significant <- (
  geo_validation$padj < 0.05 &
  abs(geo_validation$log2FoldChange) >= 1
)

geo_validation$Same_direction <- (
  sign(geo_validation$log2FoldChange) ==
  sign(geo_validation$TCGA_log2FC)
)

geo_validation$Replicated <- (
  geo_validation$GEO_significant &
  geo_validation$Same_direction
)

geo_validation <- geo_validation[
  order(geo_validation$padj),
]

geo_validation[
  geo_validation$Replicated,
]

# ============================================================
# GSE86277: TNBC specificity
# ============================================================

# 18. Download dataset
gse86277 <- getGEO(
  "GSE86277",
  GSEMatrix = TRUE
)

length(gse86277)

dim(
  exprs(gse86277[[1]])
)

# 19. Metadata
geo86277_metadata <- pData(
  gse86277[[1]]
)

table(
  geo86277_metadata$`miRNA subtype:ch1`
)

# 20. Extract IHC subtype
geo86277_metadata$IHC_subtype <- sub(
  "^ihc subtype: ",
  "",
  geo86277_metadata$characteristics_ch1.1
)

table(
  geo86277_metadata$IHC_subtype
)

# 21. Define TNBC vs non-TNBC
geo86277_metadata$TNBC_status <- ifelse(
  geo86277_metadata$IHC_subtype ==
    "Triple negative",
  "TNBC",
  "non_TNBC"
)

table(
  geo86277_metadata$TNBC_status
)

# 22. Align expression matrix
expr86277 <- exprs(
  gse86277[[1]]
)

expr86277 <- expr86277[
  ,
  match(
    geo86277_metadata$geo_accession,
    colnames(expr86277)
  )
]

all(
  colnames(expr86277) ==
    geo86277_metadata$geo_accession
)

# 23. GPL14613 annotation
platform_data <- fData(
  gse86277[[1]]
)

dim(platform_data)
colnames(platform_data)
head(platform_data, 5)

# 24. Find four replicated candidates
candidate_names <- c(
  "hsa-mir-1304",
  "hsa-mir-19b-2",
  "hsa-mir-18a",
  "hsa-mir-934"
)

for (mir in candidate_names) {

  cat("\n====================\n")
  cat(mir, "\n")

  hits <- platform_data[
    grepl(
      mir,
      platform_data$miRNA_ID_LIST,
      fixed = TRUE
    ),
    c(
      "ID",
      "miRNA_ID_LIST",
      "Sequence Type"
    )
  ]

  print(hits)
}

# 25. Candidate probe table
candidate_probe_table <- platform_data[
  platform_data$miRNA_ID_LIST %in%
    candidate_names,
  c(
    "ID",
    "miRNA_ID_LIST",
    "Sequence Type"
  )
]

candidate_probe_table

table(
  candidate_probe_table$miRNA_ID_LIST,
  candidate_probe_table$`Sequence Type`
)

# 26. Extract candidate probes
candidate_probe_ids <- c(
  "hsa-miR-1304_st",
  "hsa-miR-18a_st",
  "hsa-miR-934_st",
  "hp_hsa-mir-19b-2_st",
  "hp_hsa-mir-19b-2_x_st"
)

candidate_expr86277 <- expr86277[
  candidate_probe_ids,
  ,
  drop = FALSE
]

rownames(candidate_expr86277)

candidate_expr86277[, 1:10]

# 27. limma TNBC vs non-TNBC
library(limma)

group_86277 <- factor(
  geo86277_metadata$TNBC_status,
  levels = c(
    "non_TNBC",
    "TNBC"
  )
)

design_86277 <- model.matrix(
  ~ group_86277
)

fit_86277 <- lmFit(
  candidate_expr86277,
  design_86277
)

fit_86277 <- eBayes(
  fit_86277
)

res_86277 <- topTable(
  fit_86277,
  coef = "group_86277TNBC",
  number = Inf,
  sort.by = "P"
)

res_86277

res_86277[
  ,
  c(
    "logFC",
    "AveExpr",
    "P.Value",
    "adj.P.Val"
  )
]

# 28. Significant candidates
res_86277[
  res_86277$adj.P.Val < 0.05 &
  abs(res_86277$logFC) >= 1,
  c(
    "logFC",
    "AveExpr",
    "P.Value",
    "adj.P.Val"
  )
]

# 29. hsa-miR-934 boxplot
library(ggplot2)

miR934_plot_data <- data.frame(
  Expression = as.numeric(
    candidate_expr86277[
      "hsa-miR-934_st",
    ]
  ),
  Group = geo86277_metadata$TNBC_status
)

ggplot(
  miR934_plot_data,
  aes(
    x = Group,
    y = Expression
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.6
  ) +
  labs(
    title = "hsa-miR-934 Expression in TNBC vs Non-TNBC",
    x = "",
    y = "Expression"
  ) +
  theme_classic()

# 30. ROC analysis
library(pROC)

roc_934 <- roc(
  response = geo86277_metadata$TNBC_status,
  predictor = as.numeric(
    candidate_expr86277[
      "hsa-miR-934_st",
    ]
  ),
  levels = c(
    "non_TNBC",
    "TNBC"
  ),
  direction = "<"
)

auc_934 <- auc(roc_934)

auc_934

ci.auc(roc_934)

# 31. ROC plot
plot(
  roc_934,
  main = "ROC Curve: hsa-miR-934 for TNBC vs Non-TNBC"
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

# 32. Optimal ROC cutoff
coords_934 <- coords(
  roc_934,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity",
    "ppv",
    "npv"
  )
)

coords_934

# ============================================================
# END OF GEO ANALYSIS SO FAR
#
# Four replicated candidates:
#   hsa-mir-1304
#   hsa-mir-19b-2
#   hsa-mir-18a
#   hsa-mir-934
#
# Current lead candidate:
#   hsa-miR-934
#
# GSE86277:
#   logFC = 2.873
#   adjusted P-value = 0.0104
#   AUC = 0.8236
# ============================================================
