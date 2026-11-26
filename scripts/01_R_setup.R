library(tidyverse)
library(Seurat)
options(timeout = 600)

dir.create("results", showWarnings = FALSE)

samples  <- read_csv("data/Exp01/GSE275633_Exp01_cellbarcode_samplename.csv.gz", show_col_types = FALSE)
antigens <- read_csv("data/Exp01/GSE275633_Exp01_antigen_specificity_scores.csv.gz", show_col_types = FALSE)

table(samples$Assignment)

best <- antigens %>%
  group_by(barcode) %>%
  slice_max(antigen_specificity_score, n = 1, with_ties = FALSE) %>%
  ungroup()
view(best)
maxumi <- antigens %>%
  group_by(barcode) %>%
  summarise(max_umi = max(antigen_umi))

best <- best %>%
  left_join(samples, by = c("barcode" = "Barcode")) %>%
  left_join(maxumi, by = "barcode") %>%
  mutate(has_donor = !is.na(Assignment),
         class = case_when(
           antigen_specificity_score > 50 ~ "antigen_specific",
           max_umi < 40                   ~ "bulk",
           TRUE                           ~ "ambiguous"))

table(best$class, best$has_donor)
