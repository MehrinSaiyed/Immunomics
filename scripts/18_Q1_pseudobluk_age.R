# 18_q1_pseudobulk_age.R
# Q1: Test age-associated gene expression within broad CD8 cell states
# Input: results/ebv_annotated_harmony.rds
# Output: donor-level pseudobulk edgeR QL results

library(Seurat)
library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("results/q1_pseudobulk", showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load annotated object
# ------------------------------------------------------------

ebv <- readRDS("results/ebv_annotated_harmony.rds")

DefaultAssay(ebv) <- "RNA"

# Use original raw counts, not Harmony-corrected expression
counts <- GetAssayData(ebv, assay = "RNA", layer = "counts")

meta <- ebv@meta.data

cat("Cells loaded:", ncol(ebv), "\n")
cat("States:\n")
print(table(meta$cell_state_broad))

# ------------------------------------------------------------
# 2. Define states for primary Q1
# ------------------------------------------------------------

q1_states <- c(
  "Central_memory",
  "Effector_memory",
  "Cytotoxic_effector",
  "NK_like_cytotoxic"
)

# ------------------------------------------------------------
# 3. Build donor-by-state pseudobulk profiles
# ------------------------------------------------------------

pb_list <- list()
sample_info_list <- list()

for (state in q1_states) {
  
  cat("\nProcessing state:", state, "\n")
  
  keep_cells <- meta$cell_state_broad == state
  
  state_meta <- meta[keep_cells, , drop = FALSE]
  state_counts <- counts[, keep_cells, drop = FALSE]
  
  # Donors present in this state
  donors <- sort(unique(state_meta$donor))
  
  # Aggregate raw counts by donor
  pb <- sapply(donors, function(d) {
    rowSums(state_counts[, state_meta$donor == d, drop = FALSE])
  })
  
  # Ensure genes are rows and donors are columns
  pb <- as.matrix(pb)
  colnames(pb) <- donors
  
  # Count cells per donor in this state
  cell_n <- table(state_meta$donor)
  
  # Keep donors with at least 20 cells
  keep_donors <- names(cell_n[cell_n >= 20])
  
  cat("Donors with >=20 cells:", length(keep_donors), "\n")
  
  if (length(keep_donors) < 6) {
    cat("Skipping:", state, "- fewer than 6 usable donors\n")
    next
  }
  
  pb <- pb[, keep_donors, drop = FALSE]
  
  # Create donor metadata matching pseudobulk columns exactly
  sample_info <- state_meta %>%
    as.data.frame() %>%
    distinct(donor, age_group, experiment) %>%
    filter(donor %in% keep_donors) %>%
    arrange(match(donor, colnames(pb)))
  
  # Critical check: donor order must match exactly
  stopifnot(identical(colnames(pb), sample_info$donor))
  
  n_young <- sum(sample_info$age_group == "Y")
  n_older <- sum(sample_info$age_group == "O")
  
  cat("Young donors:", n_young, "\n")
  cat("Older donors:", n_older, "\n")
  
  # Require at least 3 donors per age group
  if (n_young < 3 | n_older < 3) {
    cat("Skipping:", state, "- insufficient donors per age group\n")
    next
  }
  
  pb_list[[state]] <- pb
  sample_info_list[[state]] <- sample_info
}

# ------------------------------------------------------------
# 4. Run edgeR QL within each state
# ------------------------------------------------------------

q1_results <- list()

for (state in names(pb_list)) {
  
  cat("\nRunning edgeR for:", state, "\n")
  
  pb <- pb_list[[state]]
  sample_info <- sample_info_list[[state]]
  
  # Ensure exact donor order
  pb <- pb[, sample_info$donor, drop = FALSE]
  
  stopifnot(identical(colnames(pb), sample_info$donor))
  
  # Create edgeR object
  y <- DGEList(
    counts = pb,
    samples = sample_info
  )
  
  # Filter low-expression genes
  keep_genes <- filterByExpr(
    y,
    group = sample_info$age_group
  )
  
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  
  # TMM normalization
  y <- normLibSizes(y)
  
  # Design: age effect adjusted for experiment
  # Reference age group is O
  # Coefficient age_groupY represents young versus older
  design <- model.matrix(
    ~ age_group + experiment,
    data = y$samples
  )
  
  cat("Design columns:\n")
  print(colnames(design))
  
  # Estimate dispersion and fit quasi-likelihood model
  y <- estimateDisp(y, design = design, robust = TRUE)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  
  # Test older versus younger
  # age_groupY is young versus older
  # Therefore, -age_groupY is older versus younger
  contrast <- makeContrasts(
    older_vs_younger = -age_groupY,
    levels = design
  )
  
  qlf <- glmQLFTest(fit, contrast = contrast)
  
  # Extract results
  results <- as.data.frame(qlf$table)
  
  # Add FDR
  results$FDR <- p.adjust(results$PValue, method = "BH")
  
  # Add gene and metadata
  results$gene <- rownames(results)
  results$state <- state
  results$n_donors <- ncol(pb)
  results$n_young_donors <- sum(sample_info$age_group == "Y")
  results$n_older_donors <- sum(sample_info$age_group == "O")
  results$comparison <- "older_vs_younger"
  
  # Ensure numeric result columns
  results$logFC <- as.numeric(results$logFC)
  results$logCPM <- as.numeric(results$logCPM)
  results$F <- as.numeric(results$F)
  results$PValue <- as.numeric(results$PValue)
  results$FDR <- as.numeric(results$FDR)
  
  q1_results[[state]] <- results
  
  write.csv(
    results,
    file.path(
      "results/q1_pseudobulk",
      paste0("q1_edgeR_", state, ".csv")
    ),
    row.names = FALSE
  )
  
  cat("Genes tested:", nrow(results), "\n")
  cat("Genes at FDR < 0.05:", sum(results$FDR < 0.05), "\n")
  cat("Genes at FDR < 0.10:", sum(results$FDR < 0.10), "\n")
}

# ------------------------------------------------------------
# 5. Combine and save all Q1 results
# ------------------------------------------------------------

q1_all <- bind_rows(q1_results)

write.csv(
  q1_all,
  "results/q1_pseudobulk/q1_edgeR_all_states.csv",
  row.names = FALSE
)

q1_significant <- q1_all %>%
  filter(FDR < 0.05)

write.csv(
  q1_significant,
  "results/q1_pseudobulk/q1_edgeR_significant_FDR0.05.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Summary table
# ------------------------------------------------------------

q1_summary <- q1_all %>%
  group_by(state) %>%
  summarise(
    genes_tested = n(),
    n_sig_FDR_0.05 = sum(FDR < 0.05),
    n_sig_FDR_0.10 = sum(FDR < 0.10),
    n_up_in_older_FDR_0.10 = sum(FDR < 0.10 & logFC > 0),
    n_down_in_older_FDR_0.10 = sum(FDR < 0.10 & logFC < 0),
    .groups = "drop"
  )

print(q1_summary)

write.csv(
  q1_summary,
  "results/q1_pseudobulk/q1_edgeR_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Top genes per state
# ------------------------------------------------------------

top_genes <- q1_all %>%
  group_by(state) %>%
  arrange(FDR, .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  top_genes,
  "results/q1_pseudobulk/q1_edgeR_top20_per_state.csv",
  row.names = FALSE
)

message("Done. Q1 pseudobulk differential expression complete.")