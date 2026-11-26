# 19_q1_stress_module_and_validation.R
# Step: validate the Q1 effector-memory stress-survival program
# Input: results/ebv_annotated_harmony.rds
# Output: stress module scores, donor-level validation, plots

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("results/q1_stress_validation", showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load object
# ------------------------------------------------------------

ebv <- readRDS("results/ebv_annotated_harmony.rds")

DefaultAssay(ebv) <- "RNA"

# ------------------------------------------------------------
# 2. Define stress-survival gene set
# ------------------------------------------------------------

stress_survival_genes <- c(
  "DDIT4",
  "DDIT3",
  "ATF3",
  "ATF4",
  "HSPA1A",
  "HSPA1B",
  "HSPH1",
  "DNAJB1",
  "BCL2",
  "BCL2L1",
  "MCL1",
  "SQSTM1",
  "BNIP3",
  "BNIP3L",
  "HIF1A",
  "VEGFA"
)

available_genes <- rownames(ebv)

stress_survival_genes <- stress_survival_genes[
  stress_survival_genes %in% available_genes
]

cat("Stress-survival genes used:\n")
print(stress_survival_genes)

# ------------------------------------------------------------
# 3. Calculate module score
# ------------------------------------------------------------

ebv <- AddModuleScore(
  object = ebv,
  features = list(stress_survival_genes),
  name = "StressSurvival"
)

meta <- ebv@meta.data

# Seurat adds the score as StressSurvival1
if (!"StressSurvival1" %in% colnames(meta)) {
  stop("StressSurvival1 was not created. Check AddModuleScore output.")
}

# ------------------------------------------------------------
# 4. Extract DDIT4 expression and metadata
# ------------------------------------------------------------

ddi_expr <- FetchData(
  ebv,
  vars = c(
    "DDIT4",
    "donor",
    "age_group",
    "experiment",
    "cell_state_broad"
  )
)

# Convert to a data frame
ddi_expr <- as.data.frame(ddi_expr)

# Check that DDIT4 was found
if (!"DDIT4" %in% colnames(ddi_expr)) {
  stop("DDIT4 was not found in the RNA assay.")
}

# ------------------------------------------------------------
# 5. Donor-level summary in effector-memory cells
# ------------------------------------------------------------

em_meta <- meta %>%
  as.data.frame() %>%
  filter(cell_state_broad == "Effector_memory")

em_stress <- em_meta %>%
  group_by(donor, age_group, experiment) %>%
  summarise(
    n_cells_stress = n(),
    mean_stress_score = mean(StressSurvival1),
    median_stress_score = median(StressSurvival1),
    .groups = "drop"
  )

em_ddit4 <- ddi_expr %>%
  filter(cell_state_broad == "Effector_memory") %>%
  group_by(donor, age_group, experiment) %>%
  summarise(
    n_cells_ddit4 = n(),
    mean_DGIT4 = mean(DDIT4),
    median_DGIT4 = median(DDIT4),
    .groups = "drop"
  )

em_donor <- em_stress %>%
  left_join(
    em_ddit4,
    by = c("donor", "age_group", "experiment")
  )

# Check that both summaries have the same donors
stopifnot(
  setequal(em_stress$donor, em_ddit4$donor)
)

write.csv(
  em_donor,
  "results/q1_stress_validation/effector_memory_donor_stress_DDIT4.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Statistical tests
# ------------------------------------------------------------

young <- em_donor$age_group == "Y"
older <- em_donor$age_group == "O"

stress_test <- wilcox.test(
  em_donor$mean_stress_score[older],
  em_donor$mean_stress_score[young]
)

ddit4_test <- wilcox.test(
  em_donor$mean_DGIT4[older],
  em_donor$mean_DGIT4[young]
)

stress_results <- data.frame(
  measure = c(
    "stress_survival_module",
    "DDIT4_expression"
  ),
  median_young = c(
    median(em_donor$mean_stress_score[young]),
    median(em_donor$mean_DGIT4[young])
  ),
  median_older = c(
    median(em_donor$mean_stress_score[older]),
    median(em_donor$mean_DGIT4[older])
  ),
  p_value = c(
    stress_test$p.value,
    ddit4_test$p.value
  )
)

stress_results$padj_BH <- p.adjust(
  stress_results$p_value,
  method = "BH"
)

print(stress_results)

write.csv(
  stress_results,
  "results/q1_stress_validation/effector_memory_stress_tests.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Plot donor-level stress score
# ------------------------------------------------------------

p1 <- ggplot(
  em_donor,
  aes(x = age_group, y = mean_stress_score, fill = age_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(color = donor),
    width = 0.15,
    size = 3,
    alpha = 0.9
  ) +
  theme_classic() +
  labs(
    x = "Age group",
    y = "Mean stress-survival module score",
    color = "Donor",
    title = "Effector-memory EBV-specific CD8 cells"
  )

ggsave(
  "results/q1_stress_validation/effector_memory_stress_score_by_donor.png",
  p1,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Plot donor-level DDIT4 expression
# ------------------------------------------------------------

p2 <- ggplot(
  em_donor,
  aes(x = age_group, y = mean_DGIT4, fill = age_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(color = donor),
    width = 0.15,
    size = 3,
    alpha = 0.9
  ) +
  theme_classic() +
  labs(
    x = "Age group",
    y = "Mean DDIT4 expression",
    color = "Donor",
    title = "Effector-memory EBV-specific CD8 cells"
  )

ggsave(
  "results/q1_stress_validation/effector_memory_DGIT4_by_donor.png",
  p2,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Save object with module score
# ------------------------------------------------------------

saveRDS(
  ebv,
  "results/ebv_annotated_harmony_stress.rds"
)

message("Done. Stress-survival validation complete.")