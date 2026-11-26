library(Seurat)
library(tidyverse)

# ---- edit ----
obj_path       <- "results/ebv_states_tcr.rds"
donor_col      <- "donor"
age_col        <- "age_group"
exp_col        <- "experiment"
ag_col         <- "antigen"
lmp2           <- "EBV_LMP2"
others         <- c("EBV_EBNA1", "EBV_EBNA3C")
survival_genes <- c(
  "BCL2", "BCL2L1", "MCL1", "DDIT3", "DDIT4", "ATF4",
  "BNIP3", "BNIP3L", "HIF1A", "VEGFA", "SQSTM1"
)
min_cells      <- 20
# --------------

obj <- readRDS(obj_path)
obj <- AddModuleScore(obj, list(survival_genes), name = "Surv", seed = 1)

g <- FetchData(obj, c("BCL2", "DDIT4"))

d <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    donor = .data[[donor_col]],
    age = .data[[age_col]],
    exp = .data[[exp_col]],
    antigen = .data[[ag_col]],
    score = Surv1
  ) %>%
  bind_cols(g) %>%
  filter(antigen %in% c(lmp2, others))

# 1. LMP2 cells per donor and experiment
cat("\n--- 1. LMP2 cells per donor and experiment ---\n")
lmp2_counts <- d %>%
  filter(antigen == lmp2) %>%
  count(age, donor, exp) %>%
  arrange(age, desc(n)) %>%
  as_tibble()
print(lmp2_counts, n = 30)

# 2. Donor-level scores per antigen (donors with >= 20 cells)
dn <- d %>%
  group_by(antigen, donor, age) %>%
  filter(n() >= min_cells) %>%
  summarise(
    score = mean(score),
    BCL2 = mean(BCL2),
    DDIT4 = mean(DDIT4),
    n = n(),
    .groups = "drop"
  )

cat("\n--- 2. Donor-level summary across antigens ---\n")
summ <- dn %>%
  group_by(antigen) %>%
  summarise(
    nY = sum(age == "Y"),
    nO = sum(age == "O"),
    median_Y = round(median(score[age == "Y"]), 4),
    median_O = round(median(score[age == "O"]), 4),
    p = if (min(nY, nO) >= 3) round(wilcox.test(score[age == "O"], score[age == "Y"])$p.value, 4) else NA_real_,
    .groups = "drop"
  )
print(summ)

# 3. LMP2 with the largest donor removed (O5 check)
top <- d %>%
  filter(antigen == lmp2) %>%
  count(donor, sort = TRUE) %>%
  slice(1) %>%
  pull(donor)

cat("\nLargest LMP2 donor:", top, "\n")
print(dn %>% filter(antigen == lmp2, donor != top) %>% count(age))

# 4. Plot across LMP2, EBNA1, and EBNA3C
p <- ggplot(dn, aes(age, score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape = donor == top), width = 0.1, size = 2.5) +
  facet_wrap(~ antigen) +
  labs(
    y = "Survival score (donor mean)",
    x = "Age group",
    shape = paste("Is", top)
  ) +
  theme_bw()

dir.create("figures", showWarnings = FALSE)
ggsave("figures/lmp2_vs_others.pdf", p, width = 8, height = 4)
print(p)