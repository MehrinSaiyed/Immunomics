# 03_clean_exp01_counts.R
# Keep RNA counts only for confident antigen-specific Exp01 cells

library(tidyverse)
library(Seurat)

exp <- "Exp01"
d   <- file.path("data", exp)

# Load cleaned metadata
clean_meta <- readRDS("results/Exp01_clean_meta.rds")

# Load full RNA count matrix
counts <- Read10X(d)
gex <- counts[["Gene Expression"]]

# Find cells present in both metadata and RNA matrix
cells <- intersect(clean_meta$barcode, colnames(gex))

cat("Metadata cells:", nrow(clean_meta), "\n")
cat("RNA matrix cells:", ncol(gex), "\n")
cat("Cells found in both:", length(cells), "\n")

# Keep only matching cells
clean_meta <- clean_meta %>%
  filter(barcode %in% cells)

gex <- gex[, clean_meta$barcode]
colnames(gex) <- clean_meta$cell_id

# Save cleaned Exp01 data
saveRDS(
  list(meta = clean_meta, gex = gex),
  "results/Exp01_clean.rds"
)

cat("Saved results/Exp01_clean.rds\n")
cat("Final cells saved:", ncol(gex), "\n")