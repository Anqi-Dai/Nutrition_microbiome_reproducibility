# Figure 2 c-h (plus the F2f / E5i raw scatters) from the cached diversity models.
#
# Consumes the fits and tidy tables cached by 10_fit_diversity_models.R; no model
# is refit here. Panels:
#   F2c  conditioning-intensity intercepts
#   F2d  food-group effects forest (main effects + abx interactions)
#   F2e  macronutrient effects forest
#   F2f  raw sweets vs ln(diversity) scatter, by antibiotic exposure
#   F2g  food-group conditional-effect curves (abx exposed vs not)
#   F2h  macronutrient conditional-effect curves
#   E5i  the pre-transplant (sdrt < 0) version of F2f, same analysis subset

source(here::here("reproduce", "human", "_human_helpers.R"))

key <- food_key()
palette <- abx_palette

fg_fit    <- readRDS(cache_path("172_fit_fg_diversity.rds"))
macro_fit <- readRDS(cache_path("172_fit_macro_diversity.rds"))
results_df       <- read_csv(cache_path("172_results_df_main_fg_diversity.csv"), show_col_types = FALSE)
results_df_macro <- read_csv(cache_path("172_results_df_macro_diversity.csv"), show_col_types = FALSE)

# F2c: intercepts -------------------------------------------------------------
intercepts <- results_df %>%
  filter(str_detect(term, "intensity")) %>%
  mutate(shortname = str_to_title(str_replace(term, "intensity", "")),
         shortname = factor(shortname, levels = c("Nonablative", "Reduced", "Ablative"))) %>%
  ggplot(aes(x = estimate, y = shortname)) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.25, linewidth = 1) +
  labs(x = "ln(diversity)", y = "") +
  scale_x_continuous(breaks = seq(1.4, 2.4, 0.2)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        panel.background = element_rect(fill = "gray96", colour = "gray96", linewidth = 0.5),
        axis.title = element_text(size = 10), aspect.ratio = 1.4)
save_panel(intercepts, "F2c_intercepts.pdf", width = 90, height = 90)

# Shared forest-plot builder for F2d / F2e ------------------------------------
# Cleans MaAsLin-style term names into readable labels, flags interactions for the
# shaded bands, and marks an effect significant when its CI does not cross zero.
forest_effects <- function(results, label_dict, level_order) {
  results %>%
    filter(effect == "fixed") %>%
    mutate(
      clean_term = term %>%
        str_replace("empiricalTRUE$", "abx") %>%
        str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") %>%
        str_replace_all(label_dict) %>%
        str_replace("empiricalTRUE:", "abx * ") %>%
        str_replace_all("_", " ")) %>%
    filter(!str_detect(clean_term, "intensity")) %>%
    mutate(is_significant = (conf.low * conf.high) >= 0,
           clean_term = factor(clean_term, levels = level_order))
}

forest_plot <- function(cleaned, band_fill, band_alpha) {
  shading <- cleaned %>%
    mutate(y_numeric = as.numeric(clean_term)) %>%
    filter(str_detect(clean_term, "\\*"))
  ggplot(cleaned, aes(x = estimate, y = clean_term)) +
    geom_rect(data = shading,
              aes(ymin = y_numeric - 0.5, ymax = y_numeric + 0.5, xmin = -Inf, xmax = Inf),
              fill = band_fill, alpha = band_alpha, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "solid", color = "blue", linewidth = 0.8) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                    size = 0.25, linewidth = 1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = "Effect Status") +
    scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                     str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
    labs(x = "ln(diversity) change", y = "") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none", axis.text.y = element_markdown(), aspect.ratio = 1.3)
}

# F2d: food-group effects -----------------------------------------------------
fg_label_dict <- setNames(key$shortname, key$fg1_name)
fg_level_order <- rev(c("abx", "EN", "TPN",
  "abx * Sweets", "Sweets", "abx * Grains", "Grains", "abx * Milk", "Milk",
  "abx * Eggs", "Eggs", "abx * Legumes", "Legumes", "abx * Meats", "Meats",
  "abx * Fruits", "Fruits", "abx * Oils", "Oils", "abx * Vegetables", "Vegetables"))
plot_fg <- forest_plot(forest_effects(results_df, fg_label_dict, fg_level_order),
                       band_fill = "#FBEADC", band_alpha = 0.7)
save_panel(plot_fg, "F2d_foodgroup_effects.pdf", width = 90, height = 140)

# F2e: macronutrient effects --------------------------------------------------
macro_dict <- c("ave_Fat_g" = "Fat", "ave_Fibers_g" = "Fiber", "ave_Sugars_g" = "Sugars")
macro_level_order <- rev(c("abx", "EN", "TPN",
  "abx * Sugars", "Sugars", "abx * Fiber", "Fiber", "abx * Fat", "Fat"))
plot_macro <- forest_plot(forest_effects(results_df_macro, macro_dict, macro_level_order),
                          band_fill = "#6B8E2340", band_alpha = 0.3)
save_panel(plot_macro, "F2e_macro_effects.pdf", width = 90, height = 140)

# Conditional-effect curves (F2g / F2h) ---------------------------------------
# Pull the food/macro x antibiotic interaction surfaces, rescale the x axis back
# to grams, and draw the predicted curve + ribbon for exposed vs not-exposed.
conditional_curves <- function(fit, scale_prefix, label_map, nrow, panel_fill, strip_fill, xlab) {
  raw <- conditional_effects(fit, surface = TRUE)
  dat <- raw[str_detect(names(raw), ":empirical")] %>%
    bind_rows(.id = "grp") %>%
    mutate(grp = str_replace(grp, ":empirical", ""),
           effect1__ = effect1__ * 100,
           across(starts_with(scale_prefix), ~ .x * 100)) %>%
    label_map()
  ggplot(dat) +
    geom_smooth(aes(x = effect1__, y = estimate__, ymin = lower__, ymax = upper__,
                    fill = effect2__, color = effect2__),
                stat = "identity", alpha = 0.3, linewidth = 1.5) +
    scale_fill_manual("antibiotics", values = palette, labels = c("not exposed", "exposed")) +
    scale_colour_manual("antibiotics", values = palette, labels = c("not exposed", "exposed")) +
    facet_wrap(~ shortname, nrow = nrow, scales = "free_x") +
    ylim(0, 3.3) +
    labs(y = "Predicted log(diversity)", x = xlab) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 8),
          panel.background = element_rect(fill = panel_fill),
          aspect.ratio = 1 / 1.5,
          strip.background = element_rect(color = "white", fill = strip_fill, linewidth = 1.5),
          axis.text = element_text(size = 8), axis.title = element_text(size = axis_text_size)) +
    guides(fill = guide_legend(title = "antibiotics"), colour = guide_legend(title = "antibiotics"))
}

fg_order <- rev(c("Vegetables", "Oils", "Fruits", "Meats", "Legumes",
                  "Eggs", "Milk", "Grains", "Sweets"))
plot_fg_curves <- conditional_curves(
  fg_fit, "fg_",
  label_map = function(df) df %>%
    left_join(key %>% select(grp = fg1_name, shortname), by = "grp") %>%
    mutate(shortname = factor(shortname, levels = fg_order)),
  nrow = 3, panel_fill = "#FAEFD140", strip_fill = "#ffecdc",
  xlab = "Food group consumed (grams)")
save_panel(plot_fg_curves, "F2g_foodgroup_curves.pdf", width = 100, height = 155)

plot_macro_curves <- conditional_curves(
  macro_fit, "ave_",
  label_map = function(df) df %>% mutate(shortname = factor(recode(grp, !!!macro_dict))),
  nrow = 1, panel_fill = "#6B8E2320", strip_fill = "#6B8E2340",
  xlab = "Macronutrient consumed (grams)")
save_panel(plot_macro_curves, "F2h_macro_curves.pdf", width = 120, height = 90)

# F2f / E5i: raw sweets scatter -----------------------------------------------
# Deterministic (no model): sweets intake vs ln(diversity), faceted by antibiotic
# exposure, with a Spearman correlation. E5i is the same plot on pre-transplant
# (sdrt < 0) samples, generated here since it is the same analysis subset.
sweets_scatter <- function(meta) {
  meta %>%
    mutate(log_div = log(simpson_reciprocal),
           antibiotics = factor(if_else(empirical == "FALSE", "not exposed", "exposed"),
                                levels = c("not exposed", "exposed"))) %>%
    ggscatter(x = "fg_sweets", y = "log_div", color = "antibiotics",
              ylab = "ln(diversity)", xlab = "Sweets consumed (grams)",
              alpha = 0.35, shape = 16, size = 1,
              add = "reg.line",
              add.params = list(color = "antibiotics", fill = "antibiotics", alpha = 0.3, size = 1.5),
              conf.int = TRUE, cor.coef = TRUE,
              cor.coeff.args = list(method = "spearman", label.sep = "\n",
                                    cor.coef.name = c("rho"), p.accuracy = 0.01, r.accuracy = 0.01,
                                    label.x.npc = "middle", label.y.npc = "top", size = 3.5)) +
    scale_x_sqrt() + scale_y_sqrt() +
    scale_fill_manual("antibiotics", values = palette, labels = c("not exposed", "exposed")) +
    scale_colour_manual("antibiotics", values = palette, labels = c("not exposed", "exposed")) +
    facet_wrap(~ antibiotics, labeller = "label_both", dir = "h") +
    theme(aspect.ratio = 1 / 1.15, legend.position = "bottom",
          legend.title = element_text(size = 9, face = "bold"),
          legend.text = element_text(size = 9),
          strip.background = element_blank(), strip.text.x = element_blank(),
          axis.text = element_text(size = 9), axis.title = element_text(size = 9)) +
    guides(fill = guide_legend(title = "antibiotics"), colour = guide_legend(title = "antibiotics"))
}

META <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)
save_panel(sweets_scatter(META), "F2f_sweets_scatter.pdf", width = 90, height = 140)
save_panel(sweets_scatter(META %>% filter(sdrt < 0)), "E5i_sweets_scatter_pretransplant.pdf",
           width = 90, height = 140)

message("F2 panels written to results/.")
