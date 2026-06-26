# Figure 4 (human): F4b, the E. faecium CLR forest.
#
# Refactor of the forest half of R63_Enterococcus_asv1_outcome.Rmd. Plot-only:
# reads the cached fixed-effect table from 40 and draws the food-group main and
# antibiotic-interaction effects on E. faecium (asv_1) CLR, the abx * Sweets
# term standing out as the one credible interval clear of zero.

source(here::here("reproduce", "human", "_human_helpers.R"))

key <- food_key()
results_df <- read_csv(cache_path("R63_results_df_asv1.csv"), show_col_types = FALSE)

# Map raw coefficient names to the panel labels (fg_* -> shortname, the
# antibiotics main/interaction terms to "abx" / "abx * <group>").
replacement_dictionary <- setNames(key$shortname, key$fg1_name)

level_order <- rev(c(
  "abx", "EN", "TPN",
  "abx * Sweets", "Sweets",
  "abx * Grains", "Grains",
  "abx * Milk", "Milk",
  "abx * Eggs", "Eggs",
  "abx * Legumes", "Legumes",
  "abx * Meats", "Meats",
  "abx * Fruits", "Fruits",
  "abx * Oils", "Oils",
  "abx * Vegetables", "Vegetables"))

cleaned_effects <- results_df |>
  filter(effect == "fixed") |>
  mutate(
    clean_term = term |>
      str_replace("empiricalTRUE$", "abx") |>
      str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") |>
      str_replace_all(replacement_dictionary) |>
      str_replace("empiricalTRUE:", "abx * ") |>
      str_replace_all("_", " ")) |>
  filter(!str_detect(clean_term, "intensity")) |>
  mutate(is_significant = (conf.low * conf.high) > 0,
         clean_term = factor(clean_term, levels = level_order))

# Shade the antibiotic-interaction rows so they read against the main effects.
shading_df <- cleaned_effects |>
  mutate(y_numeric = as.numeric(clean_term)) |>
  filter(str_detect(clean_term, "\\*"))

plot_asv1 <- ggplot(cleaned_effects, aes(x = estimate, y = clean_term)) +
  geom_rect(data = shading_df,
            aes(ymin = y_numeric - 0.5, ymax = y_numeric + 0.5, xmin = -Inf, xmax = Inf),
            fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "solid", color = "blue", linewidth = 0.8) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                  size = 0.25, linewidth = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = "Effect Status") +
  scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                   str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
  labs(title = "Enterococcus faecium", x = "CLR E. faecium change", y = "") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(face = "italic", hjust = 0.5),
        axis.text.y = element_markdown(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        aspect.ratio = 1.3)

save_panel(plot_asv1, "F4b_efaecium_forest.pdf", width = 90, height = 130)
message("F4b done")
