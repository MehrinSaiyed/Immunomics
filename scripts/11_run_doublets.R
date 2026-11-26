# 11_run_doublets.R
# Step: detect doublets in results/merged_ebv_pca.rds
# Output: results/merged_ebv_pca_doublets_flagged.rds

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)

# 1. Load normalized PCA object
ebv <- readRDS("results/merged_ebv_pca.rds")
DefaultAssay(ebv) <- "RNA"

# 2. Join Seurat v5 layers
ebv <- JoinLayers(ebv, assay = "RNA")

# Check that RNA assay now has single counts/data layers
Layers(ebv[["RNA"]])

# 3. Use experiment/capture ID for doublet detection
sample_col <- "orig.ident"

stopifnot(sample_col %in% colnames(ebv@meta.data))

table(ebv@meta.data[[sample_col]])

# 4. Convert to SingleCellExperiment
ebv_sce <- as.SingleCellExperiment(ebv, assay = "RNA")

# 5. Run scDblFinder separately by experiment/capture
ebv_sce <- scDblFinder(
  ebv_sce,
  samples = ebv_sce[[sample_col]],
  dbr = NULL,
  verbose = TRUE
)

# 6. Transfer doublet results back to Seurat metadata
doublet_results <- data.frame(
  scDblFinder.class = ebv_sce$scDblFinder.class,
  scDblFinder.score = ebv_sce$scDblFinder.score,
  row.names = colnames(ebv_sce)
)

ebv <- AddMetaData(ebv, doublet_results)

# 7. Confirm
colnames(ebv@meta.data)
table(ebv$scDblFinder.class)
summary(ebv$scDblFinder.score)

# 8. Save flagged object
saveRDS(
  ebv,
  "results/merged_ebv_pca_doublets_flagged.rds"
)

message("Done. Doublets flagged and saved; no cells removed.")