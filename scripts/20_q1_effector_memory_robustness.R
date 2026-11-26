# 20_q1_effector_memory_robustness.R
# Goal: test whether the effector-memory age signal (DDIT4 etc.) is robust.
# Checks:
#   1. Drop donors with < 20 effector-memory cells (the paper's rule)
#   2. Split the old "stress-survival" gene set into two:
#        - survival/metabolic genes (possible biology)
#        - heat-shock / immediate-early genes (typical handling artefact)
#   3. Compare young vs older donors, with and without adjusting for experiment
#   4. Check whether DDIT4 moves together with handling genes across donors
#   5. Re-run the pseudobulk differential expression (edgeR) with experiment
#      as a covariate and without the low-cell donors

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(edgeR)

dir.create("results/q1_robustness", showWarnings = FALSE)

state_of_interest <- "Effector_memory"
min_cells <- 20

# ------------------------------------------------------------
# 1. Load object and keep effector-memory cells from donors with enough cells
# ------------------------------------------------------------
ebv <- readRDS("results/ebv_annotated_harmony_stress.rds")
DefaultAssay(ebv) <- "RNA"

em <- subset(ebv, cells = colnames(ebv)[ebv$cell_state_broad == state_of_interest])

cells_per_donor <- table(em$donor)
print(cells_per_donor)

keep_donors <- names(cells_per_donor)[cells_per_donor >= min_cells]
dropped     <- setdiff(names(cells_per_donor), keep_donors)
message("Donors dropped (< ", min_cells, " cells): ",
        ifelse(length(dropped) == 0, "none", paste(dropped, collapse = ", ")))

em <- subset(em, cells = colnames(em)[em$donor %in% keep_donors])

# ------------------------------------------------------------
# 2. Two separate gene modules
# ------------------------------------------------------------
survival_genes <- c("BCL2", "BCL2L1", "MCL1", "DDIT3", "DDIT4", "ATF4",
                    "BNIP3", "BNIP3L", "HIF1A", "VEGFA", "SQSTM1")
handling_genes <- c("HSPA1A", "HSPA1B", "HSPH1", "DNAJB1", "ATF3",
                    "JUN", "FOS", "FOSB", "DUSP1")

survival_genes <- intersect(survival_genes, rownames(em))
handling_genes <- intersect(handling_genes, rownames(em))

cat("Survival genes used:", paste(survival_genes, collapse = ", "), "\n")
cat("Handling genes used:", paste(handling_genes, collapse = ", "), "\n")

em <- AddModuleScore(em, features = list(survival_genes), name = "Survival", seed = 1)
em <- AddModuleScore(em, features = list(handling_genes), name = "Handling", seed = 1)

# ------------------------------------------------------------
# 3. Donor-level table
# ------------------------------------------------------------
single_genes <- intersect(c("DDIT4", "RELB", "NFKBIA", "JUN", "FOS", "DUSP1", "HSPA1A"),
                          rownames(em))
vars <- c(single_genes, "Survival1", "Handling1")

expr <- FetchData(em, vars = c(vars, "donor", "age_group", "experiment"))

donor_tab <- expr %>%
  group_by(donor, age_group, experiment) %>%
  summarise(n_cells = n(),
            across(all_of(vars), mean),
            .groups = "drop")

donor_tab$age_group  <- factor(donor_tab$age_group, levels = c("Y", "O"))
donor_tab$experiment <- factor(donor_tab$experiment)

print(as.data.frame(donor_tab))
write.csv(donor_tab, "results/q1_robustness/effector_memory_donor_table.csv", row.names = FALSE)

# ------------------------------------------------------------
# 4. Age tests: Wilcoxon, and linear model with experiment
# ------------------------------------------------------------
measures <- intersect(c("DDIT4", "Survival1", "Handling1", "RELB", "NFKBIA"), colnames(donor_tab))

age_tests <- lapply(measures, function(m) {
  d <- donor_tab
  d$value <- d[[m]]
  w   <- wilcox.test(value ~ age_group, data = d)
  fit <- lm(value ~ experiment + age_group, data = d)
  co  <- summary(fit)$coefficients["age_groupO", ]
  data.frame(
    measure               = m,
    median_young          = median(d$value[d$age_group == "Y"]),
    median_older          = median(d$value[d$age_group == "O"]),
    wilcoxon_p            = w$p.value,
    lm_older_minus_young  = unname(co["Estimate"]),
    lm_p_adj_experiment   = unname(co["Pr(>|t|)"])
  )
})
age_tests <- do.call(rbind, age_tests)
age_tests$wilcoxon_padj_BH <- p.adjust(age_tests$wilcoxon_p, method = "BH")

print(age_tests)
write.csv(age_tests, "results/q1_robustness/effector_memory_age_tests.csv", row.names = FALSE)

# ------------------------------------------------------------
# 5. Does DDIT4 track handling genes across donors?
# ------------------------------------------------------------
cor_vars <- intersect(c("DDIT4", "JUN", "FOS", "DUSP1", "HSPA1A", "Survival1", "Handling1"),
                      colnames(donor_tab))
cor_mat <- cor(as.data.frame(donor_tab)[, cor_vars], method = "spearman")
cat("\nSpearman correlations across donors (n = ", nrow(donor_tab), "):\n", sep = "")
print(round(cor_mat, 2))
write.csv(round(cor_mat, 3), "results/q1_robustness/donor_level_correlations.csv")

# ------------------------------------------------------------
# 6. Pseudobulk differential expression with experiment covariate
# ------------------------------------------------------------
pb <- AggregateExpression(em, group.by = "donor", assays = "RNA",
                          return.seurat = FALSE)$RNA

meta_pb <- as.data.frame(donor_tab)[match(colnames(pb), donor_tab$donor), ]
stopifnot(identical(as.character(meta_pb$donor), colnames(pb)))

y <- DGEList(counts = pb)
keep <- filterByExpr(y, group = meta_pb$age_group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~ experiment + age_group, data = meta_pb)
y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "age_groupO")     # positive logFC = higher in older

tt <- topTags(qlf, n = Inf)$table
tt$gene <- rownames(tt)
write.csv(tt, "results/q1_robustness/effector_memory_DE_experiment_adjusted.csv", row.names = FALSE)

cat("\nGenes tested:", nrow(tt), "\n")
cat("Genes with FDR < 0.05:", sum(tt$FDR < 0.05), "\n")
cat("Genes with FDR < 0.10:", sum(tt$FDR < 0.10), "\n")

print(head(tt[, c("gene", "logFC", "logCPM", "PValue", "FDR")], 15))

previous_hits <- c("HDGFL3", "ANTXR2", "SKAP2", "CYB561", "GBP5",
                   "ITGA4", "ALOX5AP", "DLEU1", "LMNA", "DDIT4")
cat("\nPrevious 10 hits in this analysis:\n")
print(tt[tt$gene %in% previous_hits, c("gene", "logFC", "PValue", "FDR")])

# ------------------------------------------------------------
# 7. Plot
# ------------------------------------------------------------
plot_df <- donor_tab %>%
  select(donor, age_group, experiment, all_of(intersect(c("DDIT4", "Survival1", "Handling1"),
                                                        colnames(donor_tab)))) %>%
  pivot_longer(cols = -c(donor, age_group, experiment),
               names_to = "measure", values_to = "value")

p <- ggplot(plot_df, aes(age_group, value)) +
  geom_boxplot(aes(fill = age_group), outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(shape = experiment), width = 0.12, size = 3) +
  facet_wrap(~ measure, scales = "free_y") +
  theme_classic() +
  labs(x = "Age group", y = "Donor mean",
       title = "Effector-memory cells: donors with >= 20 cells")

ggsave("results/q1_robustness/effector_memory_robustness_plot.png",
       p, width = 9, height = 4.5, dpi = 300)

message("Done. Robustness checks complete.")