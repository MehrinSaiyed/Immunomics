# 04_clean_exp02_exp03.R
# Clean Exp02 and Exp03: confident antigen-specific cells only

library(tidyverse)
library(Seurat)

options(timeout = 600)

dir.create("results", showWarnings = FALSE)

experiments <- c("Exp02", "Exp03")

for (e in experiments) {
  
  cat("\n==============================\n")
  cat("Processing:", e, "\n")
  cat("==============================\n")
  
  d   <- file.path("data", e)
  pre <- paste0("GSE275633_", e, "_")
  
  # 1. Read metadata files
  samples <- read_csv(
    file.path(d, paste0(pre, "cellbarcode_samplename.csv.gz")),
    show_col_types = FALSE
  )
  
  antigens <- read_csv(
    file.path(d, paste0(pre, "antigen_specificity_scores.csv.gz")),
    show_col_types = FALSE
  )
  
  # 2. Keep the best antigen call per cell
  best <- antigens %>%
    group_by(barcode) %>%
    slice_max(antigen_specificity_score, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # 3. Keep only confident antigen-specific cells
  clean_meta <- best %>%
    filter(antigen_specificity_score > 50) %>%
    inner_join(samples, by = c("barcode" = "Barcode")) %>%
    separate(
      Assignment,
      into = c("donor", "sample_type", "age_group", "age", "sex", "experiment"),
      sep = "_",
      remove = FALSE
    ) %>%
    mutate(
      cell_id = paste(e, barcode, sep = "_"),
      age = as.numeric(age),
      age_group_full = case_when(
        age_group == "Y" ~ "Young",
        age_group == "O" ~ "Old",
        TRUE ~ NA_character_
      )
    ) %>%
    select(
      cell_id,
      barcode,
      antigen,
      antigen_specificity_score,
      donor,
      sample_type,
      age_group,
      age_group_full,
      age,
      sex,
      experiment,
      raw_clonotype_id,
      Assignment
    )
  
  cat("Confident antigen-specific cells:", nrow(clean_meta), "\n")
  print(table(clean_meta$donor))
  
  # 4. Save metadata
  write_csv(clean_meta, file.path("results", paste0(e, "_clean_meta.csv")))
  saveRDS(clean_meta, file.path("results", paste0(e, "_clean_meta.rds")))
  
  # 5. Load RNA counts and keep only selected cells
  counts <- Read10X(d)
  gex <- counts[["Gene Expression"]]
  
  cells <- intersect(clean_meta$barcode, colnames(gex))
  
  cat("Cells found in RNA matrix:", length(cells), "\n")
  
  clean_meta <- clean_meta %>%
    filter(barcode %in% cells)
  
  gex <- gex[, clean_meta$barcode]
  colnames(gex) <- clean_meta$cell_id
  
  # 6. Save final cleaned object
  saveRDS(
    list(meta = clean_meta, gex = gex),
    file.path("results", paste0(e, "_clean.rds"))
  )
  
  cat("Finished", e, "- Final cells saved:", ncol(gex), "\n")
  
  # 7. Free memory before next experiment
  rm(counts, gex, antigens, samples, best, clean_meta, cells)
  gc()
}