# 10_doublets_clustering.R
# Goal: flag possible doublets in the all-antigen QC object using scDblFinder.
# IMPORTANT: This script only FLAGS doublets. No cells are removed.

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

set.seed(123)

# ------------------------------------------------------------
# PART 1. Load the object
# ------------------------------------------------------------
input_file <- "results/merged_all_antigens_qc_annotated.rds"
stopifnot(file.exists(input_file))

all_qc <- readRDS(input_file)
DefaultAssay(all_qc) <- "RNA"
print(all_qc)

# ------------------------------------------------------------
# PART 2. Join the three experiment layers into one
# ------------------------------------------------------------
cells_before <- colnames(all_qc)
all_qc <- SeuratObject::JoinLayers(object = all_qc, assay = "RNA")
stopifnot(identical(cells_before, colnames(all_qc)))
print(SeuratObject::Layers(all_qc[["RNA"]]))

# ------------------------------------------------------------
# PART 3. Pull out the raw counts table (genes x cells)
# ------------------------------------------------------------
counts_matrix <- SeuratObject::LayerData(
  object = all_qc,
  assay  = "RNA",
  layer  = "counts"
)
print(dim(counts_matrix))
stopifnot(ncol(counts_matrix) == ncol(all_qc))

cell_meta <- all_qc@meta.data[colnames(counts_matrix), , drop = FALSE]
stopifnot(identical(rownames(cell_meta), colnames(counts_matrix)))

# ------------------------------------------------------------
# PART 4. Convert to SingleCellExperiment
# ------------------------------------------------------------
sce <- SingleCellExperiment(
  assays  = list(counts = counts_matrix),
  colData = S4Vectors::DataFrame(cell_meta)
)
print(sce)

# ------------------------------------------------------------
# PART 5. Check the experiment grouping
# ------------------------------------------------------------
print(
  table(
    experiment        = sce$experiment,
    virus_specificity = sce$virus_specificity,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# PART 6. Run scDblFinder
# ------------------------------------------------------------
sce <- scDblFinder(sce, samples = "experiment")

print(table(sce$scDblFinder.class))
print(table(sce$experiment, sce$scDblFinder.class))

if ("donor" %in% colnames(colData(sce))) {
  print(table(sce$donor, sce$scDblFinder.class))
}

# ------------------------------------------------------------
# PART 7. Attach calls to Seurat object & save CSV
# ------------------------------------------------------------
all_qc$scDblFinder.class <- sce$scDblFinder.class[colnames(all_qc)]
all_qc$scDblFinder.score <- sce$scDblFinder.score[colnames(all_qc)]

stopifnot(
  !anyNA(all_qc$scDblFinder.class),
  !anyNA(all_qc$scDblFinder.score)
)

doublet_calls <- data.frame(
  cell  = colnames(all_qc),
  class = as.character(all_qc$scDblFinder.class),
  score = as.numeric(all_qc$scDblFinder.score),
  stringsAsFactors = FALSE
)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(doublet_calls, "results/doublet_calls.csv", row.names = FALSE)

print(table(all_qc$scDblFinder.class, useNA = "ifany"))
message("Done. Doublets flagged and saved; no cells removed.")
