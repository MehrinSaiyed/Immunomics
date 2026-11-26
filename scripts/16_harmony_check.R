# 16_harmony_check.R
# Step: test whether experiment-level batch correction improves mixing
# Input: results/ebv_annotated.rds
# Output: Harmony UMAP and comparison plots

library(Seurat)
library(harmony)
library(ggplot2)

dir.create("results/harmony_check", showWarnings = FALSE)

# 1. Load annotated object
ebv <- readRDS("results/ebv_annotated.rds")

DefaultAssay(ebv) <- "RNA"

# 2. Run Harmony on experiment
ebv <- RunHarmony(
  ebv,
  group.by.vars = "experiment",
  reduction = "pca",
  dims.use = 1:30,
  reduction.save = "harmony"
)

# 3. UMAP using Harmony embeddings
ebv <- RunUMAP(
  ebv,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap.harmony"
)

# 4. Compare original and Harmony UMAPs
p_original <- DimPlot(
  ebv,
  reduction = "umap",
  group.by = "experiment"
) +
  ggtitle("Original PCA UMAP by experiment")

p_harmony <- DimPlot(
  ebv,
  reduction = "umap.harmony",
  group.by = "experiment"
) +
  ggtitle("Harmony UMAP by experiment")

p_combined <- p_original | p_harmony

ggsave(
  "results/harmony_check/umap_original_vs_harmony.png",
  p_combined,
  width = 14,
  height = 6,
  dpi = 300
)

# 5. Compare cell states after Harmony
p_states_original <- DimPlot(
  ebv,
  reduction = "umap",
  group.by = "cell_state_broad",
  label = TRUE
) +
  ggtitle("Original UMAP: broad states")

p_states_harmony <- DimPlot(
  ebv,
  reduction = "umap.harmony",
  group.by = "cell_state_broad",
  label = TRUE
) +
  ggtitle("Harmony UMAP: broad states")

p_states_combined <- p_states_original | p_states_harmony

ggsave(
  "results/harmony_check/umap_states_original_vs_harmony.png",
  p_states_combined,
  width = 16,
  height = 7,
  dpi = 300
)

# 6. Save Harmony object separately
saveRDS(ebv, "results/ebv_annotated_harmony.rds")

message("Done. Harmony check complete.")

ebv <- readRDS("results/ebv_annotated_harmony.rds")

DimPlot(
  ebv,
  reduction = "umap.harmony",
  group.by = "cell_state_broad",
  label = TRUE
)

DimPlot(
  ebv,
  reduction = "umap.harmony",
  group.by = "donor",
  label = FALSE
)

DimPlot(
  ebv,
  reduction = "umap.harmony",
  group.by = "age_group",
  label = FALSE
)