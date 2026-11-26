# 06_filter_qc.R
# Preliminary QC filtering and retention audit

library(Seurat)
library(tidyverse)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

obj <- readRDS("results/merged_cd8_before_qc.rds")

# Check that cell names and metadata rows are aligned
stopifnot(
  identical(colnames(obj), rownames(obj[[]])),
  !anyDuplicated(colnames(obj))
)

# Create an audit table: one row per cell
qc <- obj[[]] %>%
  as.data.frame() %>%
  rownames_to_column("qc_cell_id") %>%
  mutate(
    fail_low_genes  = nFeature_RNA < 500,
    fail_high_genes = nFeature_RNA > 6000,
    fail_low_counts = nCount_RNA < 1000,
    fail_high_counts = nCount_RNA > 25000,
    fail_high_mt   = percent.mt > 10,
    pass_qc = !(
      fail_low_genes | fail_high_genes |
        fail_low_counts | fail_high_counts | fail_high_mt
    )
  )

stopifnot(!anyNA(qc$pass_qc))

# Count individual failure flags; these may overlap
failure_counts <- qc %>%
  summarise(across(starts_with("fail_"), ~ sum(.x))) %>%
  pivot_longer(
    everything(),
    names_to = "criterion",
    values_to = "cells_flagged"
  )

# Retention per donor
donor_qc <- qc %>%
  group_by(experiment, donor, age_group_full) %>%
  summarise(
    before = n(),
    after = sum(pass_qc),
    removed = sum(!pass_qc),
    percent_retained = round(100 * after / before, 2),
    .groups = "drop"
  )

# Retention per antigen
antigen_qc <- qc %>%
  group_by(antigen) %>%
  summarise(
    before = n(),
    after = sum(pass_qc),
    removed = sum(!pass_qc),
    .groups = "drop"
  )

# Filter cells using their exact cell names
keep_cells <- qc$qc_cell_id[qc$pass_qc]
obj_qc <- subset(obj, cells = keep_cells)

# Save all antigen specificities after QC
saveRDS(obj_qc, "results/merged_all_antigens_qc.rds")

# Save a separate EBV-only object
ebv_cells <- colnames(obj_qc)[
  startsWith(as.character(obj_qc$antigen), "EBV_")
]
obj_ebv <- subset(obj_qc, cells = ebv_cells)

saveRDS(obj_ebv, "results/merged_ebv_qc.rds")

# Save audit tables
write_csv(qc, "results/QC_cell_audit.csv")
write_csv(donor_qc, "results/QC_retention_by_donor.csv")
write_csv(antigen_qc, "results/QC_retention_by_antigen.csv")

# Donor x antigen table after QC, including zero combinations
donor_antigen <- qc %>%
  group_by(donor, age_group_full, antigen) %>%
  summarise(
    cells_after_qc = sum(pass_qc),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = antigen,
    values_from = cells_after_qc,
    values_fill = 0
  )

write_csv(
  donor_antigen,
  "results/donor_antigen_after_qc.csv"
)

# Save post-QC distributions
p <- VlnPlot(
  obj_qc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "experiment",
  pt.size = 0,
  ncol = 3
)

ggsave(
  "figures/QC_after_filtering.pdf",
  plot = p,
  width = 12,
  height = 4
)

cat("\nCells before QC:", ncol(obj), "\n")
cat("Cells after QC:", ncol(obj_qc), "\n")
cat("Cells removed:", ncol(obj) - ncol(obj_qc), "\n")
cat("EBV cells after QC:", ncol(obj_ebv), "\n")
cat("EBV donors after QC:", n_distinct(obj_ebv$donor), "\n")

cat("\nIndividual QC flags (may overlap):\n")
print(failure_counts)

cat("\nQC retention by donor:\n")
print(donor_qc, n = Inf)

cat("\nQC retention by antigen:\n")
print(antigen_qc, n = Inf)

cat("\nDonor x antigen after QC:\n")
print(donor_antigen, n = Inf, width = Inf)
