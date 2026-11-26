# 02_clean_exp01.R
# Clean Exp01: keep confident antigen-specific cells with donor metadata

library(tidyverse)
library(Seurat)

exp <- "Exp01"
d   <- file.path("data", exp)
pre <- paste0("GSE275633_", exp, "_")

# 1. Read sample/donor information
samples <- read_csv(
  file.path(d, paste0(pre, "cellbarcode_samplename.csv.gz")),
  show_col_types = FALSE
)

# 2. Read antigen-specificity scores
antigens <- read_csv(
  file.path(d, paste0(pre, "antigen_specificity_scores.csv.gz")),
  show_col_types = FALSE
)

# 3. Keep the best antigen call for each cell
best_antigen <- antigens %>%
  group_by(barcode) %>%
  slice_max(antigen_specificity_score, n = 1, with_ties = FALSE) %>%
  ungroup()

# 4. Keep only confident antigen calls
confident_cells <- best_antigen %>%
  filter(antigen_specificity_score > 50)

# 5. Add donor and age metadata
meta <- confident_cells %>%
  inner_join(samples, by = c("barcode" = "Barcode")) %>%
  separate(
    Assignment,
    into = c("donor", "sample_type", "age_group", "age", "sex", "experiment"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    cell_id = paste(exp, barcode, sep = "_"),
    age = as.numeric(age),
    age_group_full = case_when(
      age_group == "Y" ~ "Young",
      age_group == "O" ~ "Old",
      TRUE ~ NA_character_
    )
  )

# 6. Make the clean metadata table
clean_meta <- meta %>%
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
    Assignment
  )

# 7. Save the metadata table
write_csv(clean_meta, "results/Exp01_clean_meta.csv")
saveRDS(clean_meta, "results/Exp01_clean_meta.rds")

# 8. Print a checkpoint
cat("Exp01 cells with confident antigen calls:", nrow(clean_meta), "\n")
cat("Donors found:\n")
print(table(clean_meta$donor))

cat("\nCells per donor and antigen:\n")
print(
  clean_meta %>%
    count(donor, age_group_full, antigen) %>%
    pivot_wider(
      names_from = antigen,
      values_from = n,
      values_fill = 0
    )
)
