library(Seurat)
library(tidyverse)

# ---- edit these lines ----
obj_path   <- "results/ebv_states_tcr.rds"
state_col  <- "cell_state_broad"
donor_col  <- "donor"
age_col    <- "age_group"
min_cells  <- 3      # smallest clone used
min_clones <- 5      # smallest number of clones per donor
drop_states <- c("Stress_handling")   # technical states; check table(d$state) and edit
# --------------------------

obj <- readRDS(obj_path)
d <- obj@meta.data %>% rownames_to_column("cell") %>%
  transmute(cell, donor = .data[[donor_col]], age = .data[[age_col]],
            state = .data[[state_col]], clone_id) %>%
  filter(!is.na(clone_id), !state %in% drop_states)
table(d$state)

# Fair measure: a clone with 20 cells will always hit more states than a clone with 3.
# So draw min_cells cells from every clone, 200 times, and count distinct states.
set.seed(1)
rare <- function(s, k, B = 200) mean(replicate(B, n_distinct(sample(s, k))))

cl <- d %>% group_by(donor, age, clone_id) %>% filter(n() >= min_cells) %>%
  summarise(n_cells = n(), states_all = n_distinct(state),
            states_at_k = rare(state, min_cells), .groups = "drop")

# donor-level table: the unit for statistics
dn <- cl %>% group_by(donor, age) %>%
  summarise(clones = n(), mean_states = mean(states_at_k), .groups = "drop")
print(dn, n = 20)

dn_use <- filter(dn, clones >= min_clones)
table(dn_use$age)
wilcox.test(mean_states ~ age, data = dn_use)

p <- ggplot(dn_use, aes(age, mean_states)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = .1, size = 2) +
  labs(y = paste("Mean states per clone at", min_cells, "cells"), x = "Age group")
dir.create("figures", showWarnings = FALSE)
ggsave("figures/clone_state_spread.pdf", p, width = 4, height = 4)
p