# 13_check_cluster_composition.R
# Step: check whether clusters are driven by experiment, age, donor, or antigen
# Input: results/ebv_singlets_clustered.rds
# Output: composition tables and plots

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("results/cluster_checks", showWarnings = FALSE)
dir.create("results/cluster_checks/plots", showWarnings = FALSE)

# 1. Load clustered singlet object
ebv <- readRDS("results/ebv_singlets_clustered.rds")

res <- "RNA_snn_res.0.6"

# 2. Ensure cluster is a factor with stable ordering
ebv@meta.data[[res]] <- factor(ebv@meta.data[[res]])

# Confirm
class(ebv@meta.data[[res]])
levels(ebv@meta.data[[res]])
table(ebv@meta.data[[res]])