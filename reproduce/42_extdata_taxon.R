# Extended Fig. E7 (human taxon family): E7a and E7e.
#
# Refactor of the heatmap and marginal-effects halves of
# R63_Enterococcus_asv1_outcome.Rmd.
#   E7a  food-group effect-size heatmap across the prevalent genera, rows
#        clustered (dendrogram) by their posterior fraction-positive profile
#   E7e  marginal E. faecium (asv_1) CLR over each food group, split by
#        antibiotic exposure, from the cached asv_1 fit's conditional effects
#
# Plot-only for E7a (reads the per-genus summary cached by 40); E7e reloads the
# cached asv_1 brms fit to pull its conditional effects.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({
  library(dendextend)
  library(ggdendro)
})

key <- food_key()
replacement_dictionary <- setNames(key$shortname, key$fg1_name)

# E7a: genus heatmap + dendrogram -----------------------------------------------
level_order <- rev(c(
  "abx * Sweets", "Sweets",
  "abx * Grains", "Grains",
  "abx * Milk", "Milk",
  "abx * Eggs", "Eggs",
  "abx * Legumes", "Legumes",
  "abx * Meats", "Meats",
  "abx * Fruits", "Fruits",
  "abx * Oils", "Oils",
  "abx * Vegetables", "Vegetables"))

# Keep only the food-group terms; assign tiered significance marks from the
# nested credible intervals (*** 99% / ** 97% / * 94%), checking the strictest
# threshold first and cascading down.
post_summary <- read_csv(cache_path("R63_genus_clr_all_models_results.csv"), show_col_types = FALSE) |>
  filter(str_detect(term, "fg_")) |>
  mutate(
    clean_term = term |>
      str_replace("empiricalTRUE$", "abx") |>
      str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") |>
      str_replace_all(replacement_dictionary) |>
      str_replace("empiricalTRUE:", "abx * ") |>
      str_replace_all("_", " "),
    mark = case_when(
      Q0.5 > 0 ~ "***",
      Q99.5 < 0 ~ "***",
      Q1.5 > 0 ~ "**",
      Q98.5 < 0 ~ "**",
      Q3 > 0 ~ "*",
      Q97 < 0 ~ "*",
      .default = ""),
    clean_term = factor(clean_term, levels = level_order))

# Cluster genera (rows) on their fraction-positive profile across the food terms.
post_matrix <- post_summary |>
  filter(!str_detect(term, "intensity")) |>
  select(outcome, term, fraction_positive) |>
  spread("term", "fraction_positive") |>
  column_to_rownames("outcome")

hc <- hclust(dist(post_matrix), method = "complete")
dendro_data <- dendro_data(hc, type = "rectangle")
ordered_genera <- hc$labels[hc$order]

dendro_plot <- ggplot(segment(dendro_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0.5)) +
  scale_x_continuous(expand = c(0, 0.5)) +
  theme_dendro()

heatmap_data <- post_summary |>
  mutate(outcome = factor(outcome, levels = ordered_genera))

heatmap_plot <- heatmap_data |>
  # Blank the non-significant tiles so the credible effects read clearly.
  mutate(Estimate = if_else(!str_detect(mark, "\\*"), 0, Estimate)) |>
  ggplot(aes(x = clean_term, y = outcome, fill = Estimate)) +
  geom_tile(color = "black", width = 0.95, height = 0.95) +
  geom_text(aes(label = mark), nudge_y = -0.15, size = 1.5, color = "black") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "firebrick", midpoint = 0) +
  labs(fill = "Effect size\n(posterior median)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown(face = "italic"),
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    legend.position = "left",
    panel.background = element_rect(colour = "white", fill = "white"),
    panel.grid = element_blank())

final_figure <- plot_grid(heatmap_plot, dendro_plot, ncol = 2,
                          rel_widths = c(1, 0.1), align = "h", axis = "tb")
save_panel(final_figure, "E7a_genus_heatmap.pdf", width = 180, height = 200)
message("E7a done")

# E7e: marginal E. faecium CLR per food group, by antibiotic exposure -----------
asv1_fit <- readRDS(cache_path("R63_fit_asv1.rds"))

# Conditional effects for the food-group x antibiotics interactions. Stretch the
# per-100g model units back to grams for the axes.
raw_conditional <- conditional_effects(asv1_fit, surface = TRUE)
condi_dat <- raw_conditional[str_detect(names(raw_conditional), ":empirical")] |>
  bind_rows(.id = "grp") |>
  mutate(grp = str_replace(grp, "\\:empirical", "")) |>
  mutate(effect1__ = effect1__ * 100,
         across(starts_with("fg_"), ~ .x * 100)) |>
  left_join(key |> select(grp = fg1_name, shortname), by = "grp") |>
  mutate(shortname = factor(shortname, levels = rev(c(
    "Vegetables", "Oils", "Fruits", "Meats",
    "Legumes", "Eggs", "Milk", "Grains", "Sweets"))))

conditional_plots <- condi_dat |>
  ggplot() +
  geom_smooth(aes(x = effect1__, y = estimate__, ymin = lower__, ymax = upper__,
                  fill = effect2__, color = effect2__),
              stat = "identity", alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual("broad-spectrum antibiotics", values = abx_palette,
                    labels = c("not exposed", "exposed")) +
  scale_colour_manual("broad-spectrum antibiotics", values = abx_palette,
                      labels = c("not exposed", "exposed")) +
  facet_wrap(~ shortname, nrow = 3, scales = "free_x") +
  labs(y = "Predicted CLR(E. faecium)", x = "Food group consumed (grams)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.key = element_rect(fill = alpha("white", 0)),
        aspect.ratio = 1,
        strip.background = element_rect(color = "white", fill = "#ffecdc", linewidth = 1.5, linetype = "solid")) +
  guides(fill = guide_legend(direction = "vertical"), color = guide_legend(direction = "vertical"))

save_panel(conditional_plots, "E7e_efaecium_marginal.pdf", width = 150, height = 160)
message("E7e done")
