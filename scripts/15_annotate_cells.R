# 14_extra_markers.R
# Step: add markers to resolve unclear clusters 9, 13, and 14
# Input: results/ebv_singlets_clustered.rds
# Output: extra marker dot plot

library(Seurat)
library(ggplot2)

dir.create("results/cluster_checks/plots", showWarnings = FALSE)

ebv <- readRDS("results/ebv_singlets_clustered.rds")

res <- "RNA_snn_res.0.6"

genes <- c(
  # MAIT-like
  "KLRB1", "SLC4A10", "ZBTB16", "CEBPD", "TRAV1-2",
  "TRDC", "TRGC1",
  
  # NK-like / cytotoxic
  "KLRC1", "KLRC2", "EOMES", "TBX21", "ZEB2",
  "KLRG1", "B3GAT1", "FCGR3A", "TYROBP", "KLRD1",
  
  # Co-stimulation / memory
  "CD27", "CD28",
  
  # Confirm prior markers
  "TCF7", "CCR7", "SELL", "IL7R", "LEF1",
  "GZMK", "GZMB", "PRF1", "NKG7", "GNLY",
  "TIGIT", "TOX", "ITGAE", "CD69",
  "JUN", "FOS", "HSPA1A"
)

genes_present <- intersect(genes, rownames(ebv))

p <- DotPlot(
  ebv,
  features = genes_present,
  group.by = res,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("Extra markers for cluster annotation")

ggsave(
  "results/cluster_checks/plots/extra_marker_dotplot_res0.6.png",
  p,
  width = 16,
  height = 8,
  dpi = 300
)

message("Done. Extra marker dot plot created.")