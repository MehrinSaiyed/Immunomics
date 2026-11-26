# 17_donor_state_composition.R
# Step: test whether broad cell-state proportions differ by age group
# Input: results/ebv_annotated_harmony.rds
# Output: donor-level state proportions and propeller results

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(speckle)

dir.create("results/state_composition", showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load annotated Harmony objecomposition_summary
wilcox_results
propeller_resultsct
# ------------------------------------------------------------

ebv <- readRDS("results/ebv_annotated_harmony.rds")

# Use the original cell-state annotation, not Harmony-corrected expression
meta <- ebv@meta.data

# Confirm inputs
table(meta$cell_state_broad)
table(meta$age_group)
table(meta$donor, meta$age_group)

# ------------------------------------------------------------
# 2. Donor-by-state cell counts
# ------------------------------------------------------------

state_counts <- meta %>%
  as.data.frame() %>%
  count(donor, age_group, cell_state_broad, name = "n_cells")

state_matrix <- state_counts %>%
  pivot_wider(
    names_from = cell_state_broad,
    values_from = n_cells,
    values_fill = 0
  )

# Put donors in rows
state_matrix <- state_matrix %>%
  arrange(donor)

donor_ids <- state_matrix$donor
age_groups <- state_matrix$age_group

count_mat <- as.matrix(state_matrix[, setdiff(colnames(state_matrix), c("donor", "age_group"))])
rownames(count_mat) <- donor_ids

# Total cells per donor
donor_totals <- rowSums(count_mat)

# Proportion of each donor's cells in each state
prop_mat <- count_mat / donor_totals

# Save tables
write.csv(
  cbind(
    donor = donor_ids,
    age_group = age_groups,
    total_cells = donor_totals,
    as.data.frame(count_mat)
  ),
  "results/state_composition/donor_state_counts.csv",
  row.names = FALSE
)

write.csv(
  cbind(
    donor = donor_ids,
    age_group = age_groups,
    total_cells = donor_totals,
    as.data.frame(prop_mat)
  ),
  "results/state_composition/donor_state_proportions.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Exclude unreliable states from primary composition test
# ------------------------------------------------------------

exclude_states <- c(
  "Donor_specific_unresolved",
  "Terminal_effector_TEMRA"
)

keep_states <- setdiff(colnames(prop_mat), exclude_states)

prop_primary <- prop_mat[, keep_states, drop = FALSE]

# ------------------------------------------------------------
# 4. Descriptive donor-level comparison
# ------------------------------------------------------------

young <- age_groups == "Y"
older <- age_groups == "O"

composition_summary <- data.frame(
  state = colnames(prop_primary),
  median_young = apply(prop_primary[young, , drop = FALSE], 2, median),
  median_older = apply(prop_primary[older, , drop = FALSE], 2, median),
  mean_young = colMeans(prop_primary[young, , drop = FALSE]),
  mean_older = colMeans(prop_primary[older, , drop = FALSE])
)

composition_summary$log2FC_median <- log2(
  (composition_summary$median_older + 1e-6) /
    (composition_summary$median_young + 1e-6)
)

composition_summary <- composition_summary %>%
  arrange(desc(abs(log2FC_median)))

print(composition_summary)

write.csv(
  composition_summary,
  "results/state_composition/state_proportion_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Wilcoxon donor-level test
# ------------------------------------------------------------

wilcox_results <- do.call(
  rbind,
  lapply(colnames(prop_primary), function(state) {
    test <- wilcox.test(
      prop_primary[older, state],
      prop_primary[young, state]
    )
    data.frame(
      state = state,
      median_young = median(prop_primary[young, state]),
      median_older = median(prop_primary[older, state]),
      p_value = test$p.value
    )
  })
)

wilcox_results$padj_BH <- p.adjust(wilcox_results$p_value, method = "BH")
wilcox_results <- wilcox_results %>%
  arrange(p_value)

print(wilcox_results)

write.csv(
  wilcox_results,
  "results/state_composition/state_proportion_wilcox.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Formal propeller test
# ------------------------------------------------------------

# propeller expects a cluster vector and sample/donor vector
cluster_vector <- meta$cell_state_broad
sample_vector <- meta$donor
group_vector <- meta$age_group

# Keep only donors with enough total cells
min_cells_per_donor <- 100
usable_donors <- names(donor_totals[donor_totals >= min_cells_per_donor])

keep_cells <- meta$donor %in% usable_donors

propeller_results <- propeller(
  clusters = cluster_vector[keep_cells],
  sample = sample_vector[keep_cells],
  group = group_vector[keep_cells]
)

print(propeller_results)

write.csv(
  as.data.frame(propeller_results),
  "results/state_composition/state_proportion_propeller.csv",
  row.names = TRUE
)

# ------------------------------------------------------------
# 7. Plot donor-level proportions
# ------------------------------------------------------------

plot_data <- as.data.frame(prop_primary)
plot_data$donor <- rownames(plot_data)
plot_data$age_group <- age_groups

plot_long <- plot_data %>%
  pivot_longer(
    cols = -c(donor, age_group),
    names_to = "state",
    values_to = "proportion"
  )

p <- ggplot(
  plot_long,
  aes(x = state, y = proportion, fill = age_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(color = age_group),
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Broad cell state",
    y = "Fraction of donor's EBV cells",
    fill = "Age group",
    color = "Age group",
    title = "Donor-level cell-state composition"
  )

ggsave(
  "results/state_composition/donor_state_proportions_boxplot.png",
  p,
  width = 12,
  height = 6,
  dpi = 300
)

message("Done. Donor-level state composition analysis complete.")