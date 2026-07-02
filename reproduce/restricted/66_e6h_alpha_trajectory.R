# Extended Data E6h: fecal microbiota alpha-diversity trajectory over the peri-transplant
# window by diet-pattern cluster (RESTRICTED). Ported from the reference
# "compare alpha trajectory.R" (the R10 cluster-trajectory family).
#
# Per-sample inverse-Simpson diversity (simpson_reciprocal) over sample day relative to
# transplant (sdrt), faceted by cluster (Cluster 1 navy, Cluster 2 gold): faint
# per-patient spaghetti under a GEE-fitted population trajectory, with a dashed line at
# each cluster's fitted nadir. The cluster (modal_diet) comes from the cleaned restricted
# df_main_clinical_outcome.rds, so this is a restricted panel that skips when absent.
#
# The reference model was log(simpson_reciprocal) ~ day_exposed_cat + modal_diet + ns(sdrt)
# and faceted the 4-way cluster x antibiotic-exposure combination. Per request the
# binarized antibiotic-exposure adjustment (day_exposed_cat) is dropped: the model is
# log(simpson_reciprocal) ~ modal_diet + ns(sdrt) and the figure is the two clusters.
# With no cluster x sdrt interaction the two fitted curves share a shape and differ by a
# constant multiplicative (log-additive) offset, the modal_diet effect.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({
  library(geepack)
  library(splines)
})

df_file <- "df_main_clinical_outcome.rds"
if (!has_restricted(df_file)) {
  message("E6h skipped: restricted df_main not found (", restricted(df_file), ").")
  message("This panel needs the diet-pattern cluster (modal_diet); place the cleaned ",
          "df_main_clinical_outcome.rds in restricted_data/.")
  quit(save = "no", status = 0)
}

results_dir <- here::here("results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cluster_cols <- c("Cluster 1" = "darkslateblue", "Cluster 2" = "darkgoldenrod2")

# Per-sample alpha diversity (released) joined to the diet-pattern cluster (restricted).
# id = numeric patient id is the GEE cluster; corstr = "exch" for repeated samples.
df_alpha <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) |>
  inner_join(read_rds(restricted(df_file)) |> select(pid, modal_diet), by = "pid") |>
  filter(!is.na(modal_diet)) |>
  mutate(modal_diet = factor(modal_diet, levels = c("Cluster 1", "Cluster 2")),
         id = parse_number(pid))

message("E6h cohort: ", n_distinct(df_alpha$pid), " patients, ", nrow(df_alpha), " samples.")

# GEE trajectory, abx-exposure term dropped. sdrt modelled with a natural spline
# (knots at the 0.1/0.5/0.9 quantiles) so the curve can dip and rebound.
fit <- geeglm(log(simpson_reciprocal) ~ modal_diet +
                ns(sdrt, knots = quantile(df_alpha$sdrt, probs = c(0.1, 0.5, 0.9))),
              id = id, corstr = "exch", data = df_alpha)
cl_row <- unlist(summary(fit)$coefficients["modal_dietCluster 2", ])
message("E6h cluster effect (Cluster 2 vs 1, log inverse-Simpson): ",
        sprintf("%.2f", cl_row["Estimate"]), " (fold ",
        sprintf("%.2f", exp(cl_row["Estimate"])), "), p = ",
        formatC(cl_row["Pr(>|W|)"], format = "g", digits = 2))

# Fitted population trajectory per cluster over the observed day range.
newdata <- expand_grid(
  modal_diet = factor(c("Cluster 1", "Cluster 2"), levels = c("Cluster 1", "Cluster 2")),
  sdrt = min(df_alpha$sdrt):max(df_alpha$sdrt))
newdata$simpson_reciprocal <- exp(predict(fit, newdata = newdata))

# Dashed reference at each cluster's fitted nadir (the diversity floor).
df_nadir <- newdata |> group_by(modal_diet) |> summarise(minval = min(simpson_reciprocal),
                                                         .groups = "drop")
message("E6h fitted nadir: ",
        paste(df_nadir$modal_diet, sprintf("%.1f", df_nadir$minval), sep = " = ", collapse = ", "))

p <- ggplot(df_alpha, aes(sdrt, simpson_reciprocal, colour = modal_diet)) +
  geom_line(aes(group = pid), alpha = 0.2, linewidth = 0.3) +
  geom_line(data = newdata, linewidth = 1.3) +
  geom_hline(data = df_nadir, aes(yintercept = minval, colour = modal_diet),
             linetype = 2, linewidth = 0.8) +
  facet_grid(. ~ modal_diet) +
  scale_colour_manual(values = cluster_cols, guide = "none") +
  scale_x_continuous(breaks = c(0, 10, 30)) +
  scale_y_continuous(breaks = c(0, 10, 20, 30)) +
  coord_cartesian(xlim = c(-10, 30), ylim = c(0, 30)) +
  labs(x = "HCT Day", y = expression("fecal microbiota " * alpha * "-diversity")) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 13))

ggsave(file.path(results_dir, "E6h_alpha_trajectory_by_cluster.pdf"), p, width = 8, height = 5)
message("Wrote results/E6h_alpha_trajectory_by_cluster.pdf")
