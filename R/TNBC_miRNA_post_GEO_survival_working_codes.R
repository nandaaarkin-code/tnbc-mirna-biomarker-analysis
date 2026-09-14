# ============================================================
# TNBC miRNA BIOMARKER PROJECT
# POST-GEO ANALYSIS: TCGA SURVIVAL ANALYSIS
# ============================================================
# This file starts AFTER the GEO validation analysis.
# It contains only the working survival-analysis code used
# for miR-934 in TCGA TNBC samples.
#
# The later combined evidence-table code that produced errors
# is intentionally NOT included.
# ============================================================


# ============================================================
# 1. INSPECT TCGA OVERALL SURVIVAL VARIABLES
# ============================================================

table(
  brca_clinical$OS_event_nature2012,
  useNA = "ifany"
)

summary(
  brca_clinical$OS_time_nature2012
)

table(
  brca_clinical$Vital_Status_nature2012,
  useNA = "ifany"
)


# ============================================================
# 2. PREPARE TCGA miR-934 EXPRESSION
# ============================================================

# Normalize counts using the DESeq2 object already created
normalized_counts <- counts(
  dds_clinical,
  normalized = TRUE
)

# Extract miR-934
mir934_normalized <- normalized_counts[
  "hsa-mir-934",
]

# Number of samples
length(mir934_normalized)


# ============================================================
# 3. PREPARE TNBC SURVIVAL METADATA
# ============================================================

tnbc_metadata <- tumor_clinical_metadata[
  tumor_clinical_metadata$clinical_group == "TNBC",
]

tnbc_patients <- tnbc_metadata$patient

# Match TNBC patients to TCGA clinical survival data
clinical_tnbc_match <- match(
  tnbc_patients,
  brca_clinical$patient_id
)

# Create survival dataset
survival_tnbc <- data.frame(
  patient = tnbc_metadata$patient,
  sample = tnbc_metadata$sample,
  stringsAsFactors = FALSE
)

# Add miR-934 expression
survival_tnbc$miR934 <- as.numeric(
  mir934_normalized[
    survival_tnbc$sample
  ]
)

# Add overall survival time
survival_tnbc$OS_time <- as.numeric(
  brca_clinical$OS_time_nature2012[
    clinical_tnbc_match
  ]
)

# Add overall survival event
survival_tnbc$OS_event <- as.numeric(
  brca_clinical$OS_event_nature2012[
    clinical_tnbc_match
  ]
)

# Add vital status
survival_tnbc$vital_status <-
  brca_clinical$Vital_Status_nature2012[
    clinical_tnbc_match
  ]

# Keep complete survival cases
survival_tnbc <- survival_tnbc[
  !is.na(survival_tnbc$miR934) &
  !is.na(survival_tnbc$OS_time) &
  !is.na(survival_tnbc$OS_event),
]

# Check final dataset
dim(survival_tnbc)

head(
  survival_tnbc
)

table(
  survival_tnbc$OS_event,
  useNA = "ifany"
)

summary(
  survival_tnbc$OS_time
)


# ============================================================
# 4. DEFINE miR-934 HIGH / LOW GROUPS
# ============================================================

# Median expression cutoff
mir934_cutoff <- median(
  survival_tnbc$miR934,
  na.rm = TRUE
)

mir934_cutoff

# Create expression groups
survival_tnbc$miR934_group <- ifelse(
  survival_tnbc$miR934 >= mir934_cutoff,
  "High",
  "Low"
)

# Convert to factor
survival_tnbc$miR934_group <- factor(
  survival_tnbc$miR934_group,
  levels = c(
    "Low",
    "High"
  )
)

# Check group sizes
table(
  survival_tnbc$miR934_group
)

# Check deaths in each group
table(
  survival_tnbc$miR934_group,
  survival_tnbc$OS_event
)


# ============================================================
# 5. KAPLAN-MEIER SURVIVAL ANALYSIS
# ============================================================

library(survival)

surv_object <- Surv(
  time = survival_tnbc$OS_time,
  event = survival_tnbc$OS_event
)

km_fit_mir934 <- survfit(
  surv_object ~ miR934_group,
  data = survival_tnbc
)

summary(
  km_fit_mir934
)


# ============================================================
# 6. LOG-RANK TEST
# ============================================================

logrank_test <- survdiff(
  surv_object ~ miR934_group,
  data = survival_tnbc
)

logrank_test

# Calculate log-rank p-value
logrank_p <- 1 - pchisq(
  logrank_test$chisq,
  df = 1
)

logrank_p


# ============================================================
# 7. KAPLAN-MEIER PLOT
# ============================================================

plot(
  km_fit_mir934,
  col = c(
    "blue",
    "red"
  ),
  lwd = 2,
  xlab = "Overall Survival (days)",
  ylab = "Survival Probability",
  main = "Overall Survival by miR-934 Expression in TNBC"
)

legend(
  "bottomleft",
  legend = c(
    "Low miR-934",
    "High miR-934"
  ),
  col = c(
    "blue",
    "red"
  ),
  lwd = 2
)

legend(
  "center",
  legend = paste(
    "Log-rank p =",
    round(logrank_p, 3)
  ),
  bty = "n"
)


# ============================================================
# 8. COX PROPORTIONAL HAZARDS MODEL
# ============================================================

cox_mir934 <- coxph(
  surv_object ~ miR934_group,
  data = survival_tnbc
)

summary(
  cox_mir934
)


# ============================================================
# 9. EXTRACT HAZARD RATIO, 95% CI AND P-VALUE
# ============================================================

cox_summary <- summary(
  cox_mir934
)

cox_results <- data.frame(
  Hazard_Ratio =
    cox_summary$coefficients[
      "miR934_groupHigh",
      "exp(coef)"
    ],
  CI_lower =
    cox_summary$conf.int[
      "miR934_groupHigh",
      "lower .95"
    ],
  CI_upper =
    cox_summary$conf.int[
      "miR934_groupHigh",
      "upper .95"
    ],
  P_Value =
    cox_summary$coefficients[
      "miR934_groupHigh",
      "Pr(>|z|)"
    ]
)

cox_results


# ============================================================
# 10. PROPORTIONAL HAZARDS ASSUMPTION
# ============================================================

ph_test <- cox.zph(
  cox_mir934
)

ph_test


# ============================================================
# 11. SURVIVAL SUMMARY TABLE
# ============================================================

survival_results_mir934 <- data.frame(
  Biomarker = "miR-934",
  Analysis = c(
    "Kaplan-Meier / Log-rank",
    "Cox proportional hazards"
  ),
  N = nrow(survival_tnbc),
  Events = sum(
    survival_tnbc$OS_event
  ),
  HR = c(
    NA,
    cox_results$Hazard_Ratio
  ),
  CI_lower = c(
    NA,
    cox_results$CI_lower
  ),
  CI_upper = c(
    NA,
    cox_results$CI_upper
  ),
  P_value = c(
    logrank_p,
    cox_results$P_Value
  ),
  stringsAsFactors = FALSE
)

survival_results_mir934


# ============================================================
# CONFIRMED RESULTS FROM THE ANALYSIS
# ============================================================
#
# TCGA TNBC survival dataset:
#   N = 135
#   Events = 24
#
# miR-934 groups:
#   Low  = 67
#   High = 68
#
# Median miR-934 cutoff:
#   39.18257
#
# Kaplan-Meier / Log-rank:
#   p = 0.8217982
#
# Cox proportional hazards:
#   HR = 0.9122953
#   95% CI = 0.4011663 - 2.074658
#   p = 0.8266687
#
# Proportional hazards assumption:
#   p = 0.74
#
# Interpretation:
#   In this TCGA TNBC cohort, miR-934 expression was NOT
#   significantly associated with overall survival.
#
# This does NOT invalidate miR-934 as a diagnostic/
# classification biomarker. Its strongest evidence in this
# project is the TNBC vs non-TNBC GEO discrimination.
#
# ============================================================
# END OF POST-GEO SURVIVAL ANALYSIS
# ============================================================
