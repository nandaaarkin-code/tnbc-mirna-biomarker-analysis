# ============================================================
# TNBC miRNA Project
# Final Candidate Comparison -> miR-934 Target Analysis
# ============================================================

# ------------------------------------------------------------
# 1. Final candidate comparison in TCGA
# ------------------------------------------------------------

final_candidates <- c(
  "hsa-mir-1304",
  "hsa-mir-19b-2",
  "hsa-mir-18a",
  "hsa-mir-934"
)

tcga_normal_final <- res_TNBC_vs_Normal[
  final_candidates,
  c("log2FoldChange", "pvalue", "padj")
]

tcga_hr_final <- res_TNBC_vs_HR[
  final_candidates,
  c("log2FoldChange", "pvalue", "padj")
]

tcga_her2_final <- res_TNBC_vs_HER2[
  final_candidates,
  c("log2FoldChange", "pvalue", "padj")
]

tcga_normal_final
tcga_hr_final
tcga_her2_final


# ------------------------------------------------------------
# 2. Extract validated miR-934 targets
# ------------------------------------------------------------

mir934_targets <- validated_df[
  validated_df$mature_mirna_id == "hsa-miR-934",
]

dim(mir934_targets)
head(mir934_targets)


# ------------------------------------------------------------
# 3. Get unique miR-934 target genes
# ------------------------------------------------------------

mir934_gene_list <- unique(
  na.omit(mir934_targets$target_symbol)
)

length(mir934_gene_list)
head(mir934_gene_list, 20)


# ------------------------------------------------------------
# 4. Convert target genes to Entrez IDs
# ------------------------------------------------------------

mir934_entrez <- bitr(
  mir934_gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

dim(mir934_entrez)
head(mir934_entrez)


# ------------------------------------------------------------
# 5. GO Biological Process enrichment
# ------------------------------------------------------------

ego_mir934 <- enrichGO(
  gene = mir934_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

as.data.frame(ego_mir934)


# ------------------------------------------------------------
# 6. GO enrichment dot plot
# ------------------------------------------------------------

dotplot(
  ego_mir934,
  showCategory = 3
)


# ------------------------------------------------------------
# 7. KEGG enrichment
# ------------------------------------------------------------

ekegg_mir934 <- enrichKEGG(
  gene = mir934_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

as.data.frame(ekegg_mir934)


# ------------------------------------------------------------
# 8. Rank miR-934 targets by validation evidence
# ------------------------------------------------------------

target_support <- as.data.frame(
  table(mir934_targets$target_symbol)
)

colnames(target_support) <- c(
  "target_symbol",
  "evidence_count"
)

target_support <- target_support[
  order(-target_support$evidence_count),
]

target_support


# ------------------------------------------------------------
# 9. Select targets with >=2 validated interaction records
# ------------------------------------------------------------

network_targets <- target_support[
  target_support$evidence_count >= 2,
]

network_targets

# Expected result:
# 21 target genes selected for the miR-934 target network
