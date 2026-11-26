library(Seurat)
library(tidyverse)

# ---- edit ----
obj_path    <- "results/ebv_states_tcr.rds"
state_col   <- "cell_state_broad"
donor_col   <- "donor"
age_col     <- "age_group"
drop_states <- c("Stress_handling")
# --------------

obj <- readRDS(obj_path)

d <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    donor = .data[[donor_col]],
    age = .data[[age_col]],
    state = .data[[state_col]],
    clone_id
  ) %>%
  filter(!is.na(clone_id), !state %in% drop_states)

# A. Does the Q3 result depend on the clone-size cutoff?
set.seed(1)
rare <- function(s, k, B = 100) mean(replicate(B, n_distinct(sample(s, k))))

spread <- function(k) {
  dn <- d %>%
    group_by(donor, age, clone_id) %>%
    filter(n() >= k) %>%
    summarise(st = rare(state, k), .groups = "drop") %>%
    group_by(donor, age) %>%
    summarise(clones = n(), m = mean(st), .groups = "drop") %>%
    filter(clones >= 5)
  
  tibble(
    k = k,
    nY = sum(dn$age == "Y"),
    nO = sum(dn$age == "O"),
    median_Y = round(median(dn$m[dn$age == "Y"]), 4),
    median_O = round(median(dn$m[dn$age == "O"]), 4),
    p = if (min(nY, nO) >= 3) round(wilcox.test(m ~ age, data = dn)$p.value, 4) else NA
  )
}

cat("\n--- PART A: Clone-size cutoff sensitivity (k = 2 to 5) ---\n")
print(map_dfr(2:5, spread))

# B. Clonality at equal cell numbers (Gini & Inverse Simpson)
gini <- function(x) {
  x <- sort(as.numeric(x))
  n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}
invsimp <- function(x) 1 / sum((x / sum(x))^2)

n_per <- d %>% count(donor)
print(n_per)

N <- min(n_per$n)
cat("\nDownsampling every donor to N =", N, "cells\n")

set.seed(1)
clon <- d %>%
  group_by(donor, age) %>%
  summarise(
    gini = mean(replicate(100, gini(table(sample(clone_id, N))))),
    invsimp = mean(replicate(100, invsimp(table(sample(clone_id, N))))),
    .groups = "drop"
  )

cat("\n--- PART B: Donor-level clonality (equal depth) ---\n")
print(clon, n = 20)

print(wilcox.test(gini ~ age, data = clon))
print(wilcox.test(invsimp ~ age, data = clon))