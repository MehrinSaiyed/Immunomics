# 12_remove_doublets_clustering.R
# Step: remove scDblFinder doublets and cluster EBV singlets
# Input: results/merged_ebv_pca_doublets_flagged.rds
# Output: results/ebv_singlets_clustered.rds

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("results/clustering", showWarnings = FALSE)
dir.create("results/clustering/plots", showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load doublet-flagged object
# ------------------------------------------------------------

ebv <- readRDS("results/merged_ebv_pca_doublets_flagged.rds")

DefaultAssay(ebv) <- "RNA"

# Confirm doublet calls
table(ebv$scDblFinder.class)

# ------------------------------------------------------------
# 2. Keep singlets only
# ------------------------------------------------------------

ebv_singlets <- subset(
  ebv,
  subset = scDblFinder.class == "singlet"
)

cat(
  "Cells before doublet removal:", ncol(ebv), "\n",
  "Cells after doublet removal:", ncol(ebv_singlets), "\n",
  "Doublets removed:", ncol(ebv) - ncol(ebv_singlets), "\n"
)

table(ebv_singlets$orig.ident)

# Save singlet object
saveRDS(ebv_singlets, "results/ebv_singlets.rds")

# ------------------------------------------------------------
# 3. Re-run normalization and PCA on singlets
# ------------------------------------------------------------

DefaultAssay(ebv_singlets) <- "RNA"

ebv_singlets <- NormalizeData(
  ebv_singlets,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

ebv_singlets <- FindVariableFeatures(
  ebv_singlets,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

ebv_singlets <- ScaleData(ebv_singlets, verbose = FALSE)
ebv_singlets <- RunPCA(ebv_singlets, npcs = 50, verbose = FALSE)

# Inspect elbow plot
p_elbow <- ElbowPlot(ebv_singlets, ndims = 50)

ggsave(
  "results/clustering/plots/elbow_plot_after_doublet_removal.png",
  p_elbow,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 4. UMAP and clustering
# ------------------------------------------------------------

# Start with 30 PCs; adjust later if the elbow plot suggests otherwise
ebv_singlets <- RunUMAP(
  ebv_singlets,
  dims = 1:30,
  reduction = "pca",
  reduction.name = "umap"
)

ebv_singlets <- FindNeighbors(
  ebv_singlets,
  dims = 1:30,
  reduction = "pca"
)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)

for (r in resolutions) {
  ebv_singlets <- FindClusters(
    ebv_singlets,
    resolution = r,
    verbose = FALSE
  )
}

saveRDS(ebv_singlets, "results/ebv_singlets_clustered.rds")

# ------------------------------------------------------------
# 5. Cluster column names
# ------------------------------------------------------------

cluster_cols <- grep(
  "^RNA_snn_res",
  colnames(ebv_singlets@meta.data),
  value = TRUE
)

cluster_cols

# ------------------------------------------------------------
# 6. UMAPs by cluster
# ------------------------------------------------------------

for (col in cluster_cols) {
  p <- DimPlot(
    ebv_singlets,
    reduction = "umap",
    group.by = col,
    label = TRUE
  ) +
    ggtitle(col)
  
  ggsave(
    file.path("results/clustering/plots", paste0("umap_", col, ".png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 7. UMAPs by metadata
# ------------------------------------------------------------

meta_plot_cols <- intersect(
  c("orig.ident", "donor", "age_group", "age", "antigen", "experiment"),
  colnames(ebv_singlets@meta.data)
)

for (col in meta_plot_cols) {
  p <- DimPlot(
    ebv_singlets,
    reduction = "umap",
    group.by = col,
    label = FALSE
  ) +
    ggtitle(col)
  
  ggsave(
    file.path("results/clustering/plots", paste0("umap_by_", col, ".png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 8. Cluster-by-donor table
# ------------------------------------------------------------

cluster_table <- ebv_singlets@meta.data %>%
  as.data.frame() %>%
  select(
    any_of(c("orig.ident", "donor", "age_group", "age", "antigen", "experiment")),
    all_of(cluster_cols)
  ) %>%
  pivot_longer(
    cols = all_of(cluster_cols),
    names_to = "cluster_column",
    values_to = "cluster"
  ) %>%
  count(
    cluster_column,
    cluster,
    orig.ident,
    donor,
    age_group,
    antigen,
    name = "n_cells"
  )

write.csv(
  cluster_table,
  "results/clustering/cluster_by_donor_table.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Marker dot plot at resolution 0.6
# ------------------------------------------------------------

markers <- c(
  "TCF7", "CCR7", "SELL", "IL7R", "LEF1",
  "GZMK", "GZMB", "PRF1", "NKG7", "GNLY", "FGFBP2",
  "CX3CR1", "ITGAE", "CD69",
  "PDCD1", "TOX", "LAG3", "TIGIT", "HAVCR2",
  "MKI67", "TOP2A", "STMN1",
  "ISG15", "IFIT1", "IFIT3", "MX1",
  "JUN", "FOS", "HSPA1A"
)

markers_present <- intersect(markers, rownames(ebv_singlets))

p_markers <- DotPlot(
  ebv_singlets,
  features = markers_present,
  group.by = "RNA_snn_res.0.6",
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("EBV CD8 T-cell markers at resolution 0.6")

ggsave(
  "results/clustering/plots/marker_dotplot_res0.6.png",
  p_markers,
  width = 14,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Top cluster markers
# ------------------------------------------------------------

all_markers <- FindAllMarkers(
  ebv_singlets,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

write.csv(
  all_markers,
  "results/clustering/all_cluster_markers_res0.6.csv",
  row.names = FALSE
)

message("Done. Doublets removed and EBV singlets clustered.")