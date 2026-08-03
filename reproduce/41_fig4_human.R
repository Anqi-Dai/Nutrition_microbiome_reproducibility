# Figure 4 (human): F4a and F4b.
#   F4a  genus-abundance vs alpha-diversity Spearman correlations, the top five
#        genera in each direction (refactor of 178_new_F4__code_for_Figure_4.Rmd)
#   F4b  the E. faecium CLR forest (refactor of the forest half of
#        R63_Enterococcus_asv1_outcome.Rmd), plot-only from the cached 40 table:
#        the abx * Sweets term stands out as the one credible interval clear of zero
#
# F4a is deterministic (Spearman); the shared correlation step lives in the helper
# (genus_diversity_spearman) so 42's E7b reuses it. F4b reads the cached fixef
# table from 40.

source(here::here("reproduce", "human", "_human_helpers.R"))

key <- food_key()

# F4a: top-5 genera in each diversity direction ---------------------------------
res <- genus_diversity_spearman()

# Keep the FDR-significant genera, take up to the five strongest correlations in
# each direction, and order them by signed rho (most positive at top after
# coord_flip). "Up to": the split is by the sign of rho, and only four genera
# correlate negatively at FDR < 0.05, so the panel draws four blue bars and five
# yellow rather than five and five.
main <- res %>%
  filter(sig05 == "FDR < 0.05") %>%
  split(.$Correlation) %>%
  map(function(df) df %>% mutate(absrho = abs(rho)) %>% slice_max(order_by = absrho, n = 5)) %>%
  bind_rows() %>%
  arrange(desc(rho))

# Styling matched to the assembled Fig. 4 artwork, which gives panel (a) a 40 x 100 pt
# bar panel on a ~10 pt genus row with the bars nearly filling each row. `hairline`
# (0.5 pt) and `genus_text_pt` (7 pt) are the shared taxon-panel constants from the
# helper, so this panel and E7a/E7b carry the same line weights and label size.
spearman_ratio <- 100 / 40

correbar <- main %>%
  mutate(genus = factor(genus, levels = main$genus)) %>%
  ggplot(aes(x = genus, y = 0, xend = genus, yend = rho, color = Correlation)) +
  # bar thickness just under the row pitch, as in the artwork
  geom_segment(linewidth = 9.5 / .pt) +
  labs(x = "", y = "Spearman correlation") +
  # three breaks, and "0" rather than "0.0": the panel is 40 pt wide
  scale_y_continuous(breaks = c(-0.2, 0, 0.4), labels = c("-0.2", "0", "0.4")) +
  scale_color_jco() +
  coord_flip() +
  theme_classic(base_size = 10) +
  theme(aspect.ratio = spearman_ratio,
        axis.text = element_text(size = genus_text_pt),
        # 4.5 pt tick labels, matching the artwork. Three of them have to fit a
        # 41 pt panel where 0.2 rho units is only 7 pt, so anything larger runs
        # the -0.2 and 0 labels together.
        axis.text.x = element_text(size = 4.5),
        axis.title = element_text(size = genus_text_pt),
        legend.position = "none",
        axis.line = element_line(linewidth = hairline),
        axis.ticks = element_line(linewidth = hairline),
        axis.text.y = element_markdown(face = "italic", size = genus_text_pt))

save_panel(correbar, "F4a_genus_diversity_spearman.pdf", width = 58, height = 44)
message("F4a done")
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
