library(Seurat)
library(tidyverse)

# 1. Load your latest annotated EBV object
obj_path  <- "results/ebv_annotated_harmony_stress.rds"
state_col <- "cell_state_broad"
donor_col <- "donor"
exp_col   <- "experiment"

obj <- readRDS(obj_path)
md  <- obj@meta.data %>% rownames_to_column("cell")
stopifnot(all(c(state_col, donor_col, exp_col, "barcode") %in% colnames(md)))

# Use the existing clean barcode column (and strip any prefix just in case)
md$bc <- sub("^.*_([ACGT]+-[0-9]+)$", "\\1", md$barcode)
print(head(md[, c("cell", "barcode", "bc")]))

# 2. Read filtered contig annotations for Exp01, Exp02, Exp03
read_tcr <- function(e) {
  f <- file.path("data", e, paste0("GSE275633_", e, "_filtered_contig_annotations.csv.gz"))
  df <- read_csv(f, show_col_types = FALSE)
  
  # First column holds the cell barcode ("Exp01_geno_barcode" in Exp01, "barcode" in others)
  raw_bc <- df[[1]]
  
  df %>%
    mutate(
      # Extract the pure 10x barcode (e.g., AAACCTG...-1) regardless of prefix
      clean_bc = sub("^.*_([ACGT]+-[0-9]+)$", "\\1", as.character(raw_bc))
    ) %>%
    filter(
      tolower(as.character(productive)) == "true",
      chain %in% c("TRA", "TRB"),
      !is.na(cdr3),
      cdr3 != ""
    ) %>%
    group_by(clean_bc, chain) %>%
    slice_max(umis, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(clean_bc, chain, cdr3, v_gene, j_gene) %>%
    pivot_wider(
      names_from = chain,
      values_from = c(cdr3, v_gene, j_gene)
    ) %>%
    mutate(!!exp_col := e)
}

tcr <- map_dfr(c("Exp01", "Exp02", "Exp03"), read_tcr)

# 3. Join TCR chains & build donor-specific clone IDs
md2 <- md %>% left_join(tcr, by = c("bc" = "clean_bc", exp_col))
stopifnot(nrow(md2) == nrow(md))

md2 <- md2 %>%
  mutate(
    has_both = !is.na(cdr3_TRA) & !is.na(cdr3_TRB),
    clone_id = ifelse(
      has_both,
      paste(.data[[donor_col]], cdr3_TRA, cdr3_TRB, sep = "|"),
      NA_character_
    )
  )

# 4. Diagnostics
cat("\n--- CHECK 1: TCR recovery per experiment ---\n")
check1 <- md2 %>%
  group_by(.data[[exp_col]]) %>%
  summarise(
    cells = n(),
    any_TRB = sum(!is.na(cdr3_TRB)),
    both_chains = sum(has_both),
    pct_both = round(100 * both_chains / cells, 1),
    .groups = "drop"
  )
print(check1)

cat("\n--- CHECK 2: Clone sizes per donor ---\n")
check2 <- md2 %>%
  filter(has_both) %>%
  count(.data[[donor_col]], clone_id) %>%
  group_by(.data[[donor_col]]) %>%
  summarise(
    clones = n(),
    cells = sum(n),
    largest = max(n),
    clones_ge3 = sum(n >= 3),
    .groups = "drop"
  )
print(check2)

# 5. Save updated object and clone table
obj$clone_id <- md2$clone_id[match(colnames(obj), md2$cell)]
saveRDS(obj, "results/ebv_states_tcr.rds")

write.csv(
  md2 %>% select(cell, all_of(c(donor_col, exp_col, state_col)), clone_id, cdr3_TRA, cdr3_TRB),
  "results/ebv_clone_table.csv",
  row.names = FALSE
)