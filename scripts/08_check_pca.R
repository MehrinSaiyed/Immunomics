# 08_check_pca.R
# Inspect PCA before choosing a clustering strategy

library(Seurat)
library(ggplot2)

obj <- readRDS("results/merged_ebv_pca.rds")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# 1. Print and save genes contributing to PC1-PC5
pc_summary <- capture.output(
  print(obj[["pca"]], dims = 1:5, nfeatures = 10)
)

cat(paste(pc_summary, collapse = "\n"), "\n")

writeLines(
  pc_summary,
  "results/PCA_top_genes_PC1_PC5.txt"
)

# 2. Confirm antigen composition
cat("\nAntigen counts by experiment:\n")
print(table(obj$antigen, obj$experiment))

# 3. Plot PCA by donor, experiment, and age
groups <- c("donor", "experiment", "age_group_full")

pc_pairs <- list(
  PC1_PC2 = c(1, 2),
  PC3_PC4 = c(3, 4),
  PC4_PC5 = c(4, 5)
)

for (group in groups) {
  for (pair_name in names(pc_pairs)) {
    
    p <- DimPlot(
      obj,
      reduction = "pca",
      dims = pc_pairs[[pair_name]],
      group.by = group,
      shuffle = TRUE,
      seed = 123,
      pt.size = 0.2
    )
    
    ggsave(
      filename = paste0(
        "figures/PCA_", pair_name, "_by_", group, ".pdf"
      ),
      plot = p,
      width = 9,
      height = 6
    )
  }
}

# 4. Explore selected differentiation-associated genes
marker_genes <- intersect(
  c("TCF7", "LEF1", "CCR7", "GZMB", "NKG7", "KLRD1"),
  rownames(obj)
)

p_markers <- FeaturePlot(
  obj,
  reduction = "pca",
  features = marker_genes,
  ncol = 3,
  order = TRUE
)

ggsave(
  "figures/PCA_marker_expression.pdf",
  plot = p_markers,
  width = 12,
  height = 8
)

# 5. Inspect relationships between PCs and QC metrics
pc_scores <- Embeddings(obj, reduction = "pca")[, 1:10, drop = FALSE]

qc_metrics <- obj[[]][
  rownames(pc_scores),
  c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  drop = FALSE
]

pc_qc <- cor(
  pc_scores,
  as.matrix(qc_metrics),
  method = "spearman",
  use = "complete.obs"
)

cat("\nPC correlations with QC metrics:\n")
print(round(pc_qc, 3))

write.csv(
  pc_qc,
  "results/PCA_QC_correlations.csv"
)

cat("\nPCA checks finished. No clustering performed yet.\n")