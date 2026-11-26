library(Seurat)
library(tidyverse)

obj_path   <- "results/ebv_states_tcr.rds"
state_col  <- "cell_state_broad"
em_label   <- "Effector_memory"
donor_col  <- "donor"
age_col    <- "age_group"
exp_col    <- "experiment"
ag_col     <- "antigen"
min_cells  <- 20

survival_genes <- c(
  "BCL2", "BCL2L1", "MCL1", "DDIT3", "DDIT4", "ATF4",
  "BNIP3", "BNIP3L", "HIF1A", "VEGFA", "SQSTM1"
)
handling_genes <- c(
  "FOS", "FOSB", "JUN", "JUNB", "HSPA1A", "HSPA1B",
  "DNAJB1", "IER2", "DUSP1", "HSPH1"
)

obj <- readRDS(obj_path)
em  <- subset(obj, cells = colnames(obj)[obj@meta.data[[state_col]] == em_label])
em  <- AddModuleScore(em, list(survival_genes), name = "Surv", seed = 1)
em  <- AddModuleScore(em, list(intersect(handling_genes, rownames(em))), name = "Hand", seed = 1)

d <- em@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    donor = .data[[donor_col]],
    age = .data[[age_col]],
    antigen = .data[[ag_col]],
    experiment = .data[[exp_col]],
    clone_id,
    score = Surv1,
    hand = Hand1
  )

test_score <- function(df, label) {
  dn <- df %>%
    group_by(donor, age) %>%
    filter(n() >= min_cells) %>%
    summarise(s = mean(score), .groups = "drop")
  
  nY <- sum(dn$age == "Y")
  nO <- sum(dn$age == "O")
  p  <- if (nY >= 3 && nO >= 3) wilcox.test(s ~ age, data = dn)$p.value else NA
  
  tibble(
    check = label,
    nY, nO,
    median_Y = round(median(dn$s[dn$age == "Y"]), 4),
    median_O = round(median(dn$s[dn$age == "O"]), 4),
    p = round(p, 4)
  )
}

set.seed(1)
res <- bind_rows(
  test_score(d, "all cells"),
  map_dfr(unique(d$donor), ~ test_score(filter(d, donor != .x), paste("drop", .x))),
  test_score(
    d %>% filter(!is.na(clone_id)) %>% group_by(donor, clone_id) %>% slice_sample(n = 1) %>% ungroup(),
    "one cell per clone"
  ),
  test_score(filter(d, hand < quantile(hand, 0.9)), "drop top 10% handling")
)

print(res, n = 40)
write.csv(res, "results/q1_robustness.csv", row.names = FALSE)

# Donor-level linear model adjusting for experiment (requires no extra packages)
dn_model <- d %>%
  group_by(donor, age, experiment) %>%
  filter(n() >= min_cells) %>%
  summarise(score = mean(score), .groups = "drop") %>%
  mutate(age = factor(age, levels = c("Y", "O")))

m_lm <- lm(score ~ age + experiment, data = dn_model)
print(summary(m_lm)$coefficients)