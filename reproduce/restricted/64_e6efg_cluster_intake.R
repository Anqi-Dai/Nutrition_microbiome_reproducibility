# Extended Data E6 e,f,g: daily dietary intake over the peri-transplant window by
# diet-pattern cluster (RESTRICTED). Ported from R10_patients_trajectory_clusters.Rmd.
#
#   E6e  calories       (Daily intake, 10^3 kcal)
#   E6f  sugars         (Daily intake, proportion of macro mass)
#   E6g  other carbs / fat / protein (Daily intake, proportion)
#
# Per-patient-day macro totals come from the released diet table; the diet-pattern
# cluster (modal_diet) comes from the cleaned restricted df_main_clinical_outcome.rds,
# so this is a restricted panel. Skips cleanly when df_main is absent.
#
# Each panel: boxplots per HCT day -7..12, faceted by cluster (Cluster 1 blue left,
# Cluster 2 orange right), red median-trajectory line, and a between-cluster GEE
# Wald p-value (proportion or kcal ~ cluster + ns(day); exchangeable, id = patient).

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({
  library(geepack)
  library(splines)
})

df_file <- "df_main_clinical_outcome.rds"
if (!has_restricted(df_file)) {
  message("E6 e,f,g skipped: restricted df_main not found (", restricted(df_file), ").")
  message("This panel needs the diet-pattern cluster (modal_diet); place the cleaned ",
          "df_main_clinical_outcome.rds in restricted_data/.")
  quit(save = "no", status = 0)
}

results_dir <- here::here("results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cluster_cols <- c("Cluster 1" = "darkslateblue", "Cluster 2" = "darkgoldenrod2")

# Diet-pattern cluster per patient (restricted).
clusters <- read_rds(restricted(df_file)) |> select(pid, modal_diet)

# Per-patient-day macro totals from the released diet table; keep the day -7..12
# window and drop zero-calorie days (spurious proportions).
dens <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE) |>
  group_by(pid, fdrt) |>
  summarise(d_carb          = sum(Carbohydrates_g),
            d_Protein_g     = sum(Protein_g),
            d_Fat_g         = sum(Fat_g),
            d_Calories_kcal = sum(Calories_kcal),
            d_Sugars_g      = sum(Sugars_g), .groups = "drop") |>
  filter(d_Calories_kcal > 0, fdrt %in% -7:12) |>
  inner_join(clusters, by = "pid") |>
  mutate(macro_sum       = d_carb + d_Protein_g + d_Fat_g,
         cal_k           = d_Calories_kcal / 1000,
         sugar_prop      = if_else(macro_sum > 0, d_Sugars_g / macro_sum, 0),
         other_carb_prop = if_else(macro_sum > 0, (d_carb - d_Sugars_g) / macro_sum, 0),
         fat_prop        = if_else(macro_sum > 0, d_Fat_g / macro_sum, 0),
         protein_prop    = if_else(macro_sum > 0, d_Protein_g / macro_sum, 0),
         modal_diet      = factor(modal_diet, levels = c("Cluster 1", "Cluster 2")))

# Between-cluster GEE Wald p on the raw (per-day) values, day modelled with a
# natural spline (knots at the 0.2/0.8 quantiles), exchangeable working correlation.
gee_p <- function(yvar) {
  d <- dens |> mutate(fdrt_n = fdrt)
  f <- reformulate(c("modal_diet", "ns(fdrt_n, knots = quantile(fdrt_n, c(0.2, 0.8)))"),
                   response = yvar)
  m <- geeglm(f, id = factor(pid), corstr = "exch", data = d)
  summary(m)$coefficients[2, 4]
}
fmt_p <- function(p) if (p < 0.001) "p<0.001" else paste0("p=", formatC(p, format = "g", digits = 2))

# A between-cluster significance bracket drawn as its own thin plot. Stacking it
# over the boxplot with align = "v" / axis = "lr" lines up the two plots' left and
# right edges, so the bracket's x in [0,1] spans the two-facet panel area and its
# ends at x = 0.25 / 0.75 land on the centres of the Cluster 1 (left) and Cluster 2
# (right) facets. (No aspect.ratio on the boxplot, so the facets fill their halves
# and those centres stay at 0.25 / 0.75.)
bracket_plot <- function(pval) {
  ggplot() +
    annotate("segment", x = 0.25, xend = 0.75, y = 1, yend = 1,   linewidth = 0.4) +
    annotate("segment", x = 0.25, xend = 0.25, y = 1, yend = 0.7, linewidth = 0.4) +
    annotate("segment", x = 0.75, xend = 0.75, y = 1, yend = 0.7, linewidth = 0.4) +
    annotate("text",    x = 0.50, y = 1.6, label = pval, size = 3) +
    scale_x_continuous(limits = c(0, 1),   expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 2.3), expand = c(0, 0)) +
    theme_void()
}

# One faceted boxplot panel (grey title strip) topped by its between-cluster bracket.
intake_panel <- function(yvar, title, ylab, ylim = NULL, pval) {
  fdf <- dens |> mutate(fdrt = factor(fdrt))
  p <- ggplot(fdf, aes(fdrt, .data[[yvar]], fill = modal_diet, colour = modal_diet)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
    stat_summary(fun = median, geom = "line", aes(group = 1), colour = "red", linewidth = 0.5) +
    scale_x_discrete(breaks = seq(0, 10, 10)) +
    scale_colour_manual(values = cluster_cols) +
    scale_fill_manual(values = cluster_cols) +
    facet_wrap(~ modal_diet) +
    labs(title = title, x = "HCT Day", y = ylab) +
    theme_bw() +
    theme(legend.position = "none", strip.text = element_blank(),
          plot.title = ggtext::element_textbox_simple(
            halign = 0.5, fill = "grey85", padding = margin(2, 0, 2, 0),
            margin = margin(b = 1), size = 9))
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  plot_grid(bracket_plot(pval), p, ncol = 1, align = "v", axis = "lr",
            rel_heights = c(0.13, 1))
}

p_cal   <- gee_p("cal_k")
p_sugar <- gee_p("sugar_prop")
p_ocarb <- gee_p("other_carb_prop")
p_fat   <- gee_p("fat_prop")
p_prot  <- gee_p("protein_prop")
message("E6 GEE p-values: calories ", fmt_p(p_cal), "; sugars ", fmt_p(p_sugar),
        "; other carbs ", fmt_p(p_ocarb), "; fat ", fmt_p(p_fat), "; protein ", fmt_p(p_prot))

# E6e calories
e6e <- intake_panel("cal_k", "calories", expression("Daily intake (" * 10^3 ~ "kcal)"),
                    pval = fmt_p(p_cal))
ggsave(file.path(results_dir, "E6e_calories_by_cluster.pdf"), e6e, width = 2.8, height = 4)

# E6f sugars
e6f <- intake_panel("sugar_prop", "sugars", "Daily intake (proportion)",
                    ylim = c(0, 1), pval = fmt_p(p_sugar))
ggsave(file.path(results_dir, "E6f_sugars_by_cluster.pdf"), e6f, width = 2.8, height = 4)

# E6g other carbs / fat / protein (each with its own between-cluster bracket)
e6g <- plot_grid(
  intake_panel("other_carb_prop", "other carbs", "Daily intake (proportion)", ylim = c(0, 1), pval = fmt_p(p_ocarb)),
  intake_panel("fat_prop",        "fat",         "Daily intake (proportion)", ylim = c(0, 1), pval = fmt_p(p_fat)),
  intake_panel("protein_prop",    "protein",     "Daily intake (proportion)", ylim = c(0, 1), pval = fmt_p(p_prot)),
  nrow = 1)
ggsave(file.path(results_dir, "E6g_macros_by_cluster.pdf"), e6g, width = 7.5, height = 4)

message("Wrote E6 e,f,g panels to results/ (", n_distinct(dens$pid), " patients).")
