# 07_normalize_pca.R
# Prepare the EBV-only QC object for exploratory analysis

library(Seurat)
library(ggplot2)

set.seed(123)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

obj <- readRDS("results/merged_ebv_qc.rds")
DefaultAssay(obj) <- "RNA"

cat("EBV cells:", ncol(obj), "\n")
cat("Donors:", length(unique(obj$donor)), "\n")

# Record software versions
writeLines(
  capture.output(sessionInfo()),
  "results/sessionInfo_preprocessing.txt"
)

# Log-normalize RNA expression
obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# Select variable genes
obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

features <- VariableFeatures(obj)

cat("Variable genes selected:", length(features), "\n")

# Scale only selected genes to limit memory use
obj <- ScaleData(
  obj,
  features = features,
  verbose = FALSE
)

# Principal component analysis
obj <- RunPCA(
  obj,
  features = features,
  npcs = 40,
  verbose = FALSE
)

# Save first so results are preserved before plotting
saveRDS(obj, "results/merged_ebv_pca.rds")

p_elbow <- ElbowPlot(obj, ndims = 40)

ggsave(
  "figures/PCA_elbow.pdf",
  plot = p_elbow,
  width = 7,
  height = 5
)

p_exp <- DimPlot(
  obj,
  reduction = "pca",
  group.by = "experiment"
)

ggsave(
  "figures/PCA_by_experiment.pdf",
  plot = p_exp,
  width = 7,
  height = 5
)

p_age <- DimPlot(
  obj,
  reduction = "pca",
  group.by = "age_group_full"
)

ggsave(
  "figures/PCA_by_age.pdf",
  plot = p_age,
  width = 7,
  height = 5
)

p_antigen <- DimPlot(
  obj,
  reduction = "pca",
  group.by = "antigen"
)

ggsave(
  "figures/PCA_by_antigen.pdf",
  plot = p_antigen,
  width = 9,
  height = 5
)

cat("\nNormalization and PCA finished.\n")
cat("Saved results/merged_ebv_pca.rds\n")