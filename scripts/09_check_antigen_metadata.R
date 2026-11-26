# 09_check_antigen_metadata.R
# Purpose:
#   Verify the existing EBV subset.
#   Add an explicit antigen-derived virus-specificity field.
#   Preserve all original files and metadata columns.

library(Seurat)

# ------------------------------------------------------------
# 1. Load existing objects
# ------------------------------------------------------------

all_qc <- readRDS("results/merged_all_antigens_qc.rds")
ebv <- readRDS("results/merged_ebv_pca.rds")

# ------------------------------------------------------------
# 2. Derive virus specificity from antigen labels
# ------------------------------------------------------------

derive_virus <- function(antigen) {
  antigen <- as.character(antigen)
  
  virus <- rep(NA_character_, length(antigen))
  
  is_ebv <- !is.na(antigen) & startsWith(antigen, "EBV_")
  is_vzv <- !is.na(antigen) & startsWith(antigen, "VZV_")
  
  virus[is_ebv] <- "EBV"
  virus[is_vzv] <- "VZV"
  
  virus
}

all_qc$virus_specificity <- derive_virus(all_qc$antigen)
ebv$virus_specificity <- derive_virus(ebv$antigen)

# ------------------------------------------------------------
# 3. Verify that the existing EBV subset matches exactly
# ------------------------------------------------------------

stopifnot(!anyNA(all_qc$virus_specificity))
stopifnot(!anyNA(ebv$virus_specificity))

expected_ebv_cells <- colnames(all_qc)[
  all_qc$virus_specificity == "EBV"
]

stopifnot(
  setequal(colnames(ebv), expected_ebv_cells),
  all(ebv$virus_specificity == "EBV")
)

# ------------------------------------------------------------
# 4. Create and save verification tables
# ------------------------------------------------------------

antigen_counts <- as.data.frame(
  table(
    virus_specificity = all_qc$virus_specificity,
    antigen = all_qc$antigen
  ),
  stringsAsFactors = FALSE
)

antigen_counts <- antigen_counts[
  antigen_counts$Freq > 0,
]

names(antigen_counts)[3] <- "n_cells"

experiment_virus_counts <- as.data.frame(
  table(
    experiment = all_qc$experiment,
    virus_specificity = all_qc$virus_specificity
  ),
  stringsAsFactors = FALSE
)

names(experiment_virus_counts)[3] <- "n_cells"

print(antigen_counts)
print(experiment_virus_counts)

write.csv(
  antigen_counts,
  "results/antigen_specificity_counts.csv",
  row.names = FALSE
)

write.csv(
  experiment_virus_counts,
  "results/experiment_virus_specificity_counts.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Save annotated copies; do not overwrite the originals
# ------------------------------------------------------------

saveRDS(
  all_qc,
  "results/merged_all_antigens_qc_annotated.rds"
)

saveRDS(
  ebv,
  "results/merged_ebv_pca_annotated.rds"
)

message(
  "Verification passed: the existing EBV subset matches ",
  "the EBV antigen labels in the broader QC object."
)
