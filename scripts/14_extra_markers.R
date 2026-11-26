# 14_extra_markers.R
# Step: add markers to resolve unclear clusters and check alternative identities
# Input: results/ebv_singlets_clustered.rds
# Output: extra marker dot plot and marker table

library(Seurat)
library(dplyr)
library(ggplot2)

dir.create("results/cluster_checks", showWarnings = FALSE)
dir.create("results/cluster_checks/plots", showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load clustered singlet object
# ------------------------------------------------------------

ebv <- readRDS("results/ebv_singlets_clustered.rds")

res <- "RNA_snn_res.0.6"

# Ensure cluster is a factor
ebv@meta.data[[res]] <- factor(ebv@meta.data[[res]])

# Confirm clusters
table(ebv@meta.data[[res]])

# ------------------------------------------------------------
# 2. Define marker sets
# ------------------------------------------------------------

# MAIT-like markers
mait_markers <- c(
  "KLRB1", "SLC4A10", "ZBTB16", "CEBPD",
  "TRAV1-2", "TRDC", "TRGC1",
  "IL18R1", "CXCR6", "RORC"
)

# NK-like / innate cytotoxic markers
nk_markers <- c(
  "KLRC1", "KLRC2", "FCGR3A", "TYROBP", "KLRD1",
  "NKG7", "GNLY", "PRF1", "GZMB"
)

# Effector / TEMRA markers
effector_markers <- c(
  "EOMES", "TBX21", "ZEB2", "KLRG1", "B3GAT1",
  "GZMB", "PRF1", "GNLY", "FGFBP2", "CX3CR1"
)

# Memory / stem-like markers
memory_markers <- c(
  "TCF7", "CCR7", "SELL", "IL7R", "LEF1",
  "CD27", "CD28"
)

# Exhaustion and tissue-resident markers
dysfunction_markers <- c(
  "PDCD1", "TOX", "LAG3", "TIGIT", "HAVCR2",
  "ITGAE", "CD69"
)

# Stress / immediate-early markers
stress_markers <- c(
  "JUN", "FOS", "HSPA1A", "HSPA1B", "FOSB",
  "DUSP1", "ATF3"
)

# LMP2-relevant markers
lmp2_markers <- c(
  "GZMK", "TIGIT", "CD27", "CD28",
  "TCF7", "IL7R", "CCR7"
)

# ------------------------------------------------------------
# 3. Keep only genes present in the object
# ------------------------------------------------------------

all_markers <- unique(
  c(
    mait_markers,
    nk_markers,
    effector_markers,
    memory_markers,
    dysfunction_markers,
    stress_markers,
    lmp2_markers
  )
)

genes_present <- intersect(all_markers, rownames(ebv))

cat(
  "Requested markers:", length(all_markers), "\n",
  "Markers present in object:", length(genes_present), "\n"
)

print(setdiff(all_markers, genes_present))

# ------------------------------------------------------------
# 4. Extra marker dot plot
# ------------------------------------------------------------

p_extra <- DotPlot(
  ebv,
  features = genes_present,
  group.by = res,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("Extra markers for cluster annotation at resolution 0.6")

ggsave(
  "results/cluster_checks/plots/extra_marker_dotplot_res0.6.png",
  p_extra,
  width = 18,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 5. Focused dot plots
# ------------------------------------------------------------

# MAIT-like check
p_mait <- DotPlot(
  ebv,
  features = intersect(mait_markers, rownames(ebv)),
  group.by = res,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("MAIT-like markers")

ggsave(
  "results/cluster_checks/plots/mait_marker_dotplot_res0.6.png",
  p_mait,
  width = 10,
  height = 6,
  dpi = 300
)

# NK-like / cytotoxic check
p_nk <- DotPlot(
  ebv,
  features = intersect(nk_markers, rownames(ebv)),
  group.by = res,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("NK-like and cytotoxic markers")

ggsave(
  "results/cluster_checks/plots/nk_cytotoxic_marker_dotplot_res0.6.png",
  p_nk,
  width = 10,
  height = 6,
  dpi = 300
)

# Effector / TEMRA check
p_effector <- DotPlot(
  ebv,
  features = intersect(effector_markers, rownames(ebv)),
  group.by = res,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("Effector and TEMRA markers")

ggsave(
  "results/cluster_checks/plots/effector_marker_dotplot_res0.6.png",
  p_effector,
  width = 10,
  height = 6,
  dpi = 300
)

# Memory / stem-like check
p_memory <- DotPlot(
  ebv,
  features = intersect(memory_markers, rownames(ebv)),
  group.by = res,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("Memory and stem-like markers")

ggsave(
  "results/cluster_checks/plots/memory_marker_dotplot_res0.6.png",
  p_memory,
  width = 10,
  height = 6,
  dpi = 300
)

# Stress check
p_stress <- DotPlot(
  ebv,
  features = intersect(stress_markers, rownames(ebv)),
  group.by = res,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("Stress and immediate-early markers")

ggsave(
  "results/cluster_checks/plots/stress_marker_dotplot_res0.6.png",
  p_stress,
  width = 10,
  height = 6,
  dpi = 300
)

# LMP2-relevant check
p_lmp2 <- DotPlot(
  ebv,
  features = intersect(lmp2_markers, rownames(ebv)),
  group.by = res,
  scale = TRUE
) +
  RotatedAxis() +
  ggtitle("LMP2-relevant markers")

ggsave(
  "results/cluster_checks/plots/lmp2_marker_dotplot_res0.6.png",
  p_lmp2,
  width = 10,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 6. Save marker expression summary
# ------------------------------------------------------------

marker_expression <- DotPlot(
  ebv,
  features = genes_present,
  group.by = res,
  scale = TRUE
)$data

write.csv(
  marker_expression,
  "results/cluster_checks/extra_marker_expression_res0.6.csv",
  row.names = FALSE
)

message("Done. Extra marker plots and expression table created.")