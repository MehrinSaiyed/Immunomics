# 05_merge_qc.R
# Merge Exp01, Exp02 and Exp03; calculate QC metrics; filter low-quality cells

library(tidyverse)
library(Seurat)

exps <- c("Exp01", "Exp02", "Exp03")

# 1. Combine metadata from all experiments
meta_all <- map_dfr(exps, function(e) {
  readRDS(file.path("results", paste0(e, "_clean.rds")))$meta
})

# 2. Check the merged metadata
cat("Total cells:", nrow(meta_all), "\n")
cat("Total donors:", n_distinct(meta_all$donor), "\n\n")

cat("Donors by age group:\n")
print(
  meta_all %>%
    distinct(donor, age_group_full, age, sex, experiment) %>%
    arrange(age_group_full, donor)
)

cat("\nCells per donor:\n")
print(table(meta_all$donor))

cat("\nCells per antigen:\n")
print(table(meta_all$antigen))

# 3. Create Seurat objects for each experiment
obj_list <- map(exps, function(e) {
  x <- readRDS(file.path("results", paste0(e, "_clean.rds")))
  
  CreateSeuratObject(
    counts = x$gex,
    meta.data = x$meta,
    project = e
  )
})

# 4. Merge into one Seurat object
obj <- merge(obj_list[[1]], y = obj_list[-1])

# 5. Calculate QC metrics
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# 6. Inspect QC before filtering
pdf("figures/QC_before_filtering.pdf", width = 12, height = 4)
print(VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

# 7. Save unfiltered merged object
saveRDS(obj, "results/merged_cd8_before_qc.rds")

# 8. Print QC summaries
cat("\nQC summary before filtering:\n")
print(summary(obj$nFeature_RNA))
print(summary(obj$nCount_RNA))
print(summary(obj$percent.mt))

cat("\nMerged object saved: results/merged_cd8_before_qc.rds\n")


# Checking quantiles

obj <- readRDS("results/merged_cd8_before_qc.rds")

quantile(obj$nFeature_RNA, c(0.01, 0.05, 0.5, 0.95, 0.99))
quantile(obj$nCount_RNA, c(0.01, 0.05, 0.5, 0.95, 0.99))
quantile(obj$percent.mt, c(0.01, 0.05, 0.5, 0.95, 0.99))

table(obj$nFeature_RNA < 500)
table(obj$percent.mt > 10)
table(obj$nCount_RNA > 25000)