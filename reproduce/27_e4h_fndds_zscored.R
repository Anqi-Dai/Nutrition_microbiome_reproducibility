# Extended Fig. E4h: the F2d food-group diversity forest, but with the FNDDS food
# groups Z-SCORED (standardised to mean 0, SD 1) rather than expressed per 100 g, so
# every food-group effect reads as ln(diversity) change per 1 SD of intake.
#
# Same model and predictor set as the F2d food-group diversity model
# (10_fit_diversity_models.R / 11_fig2_diversity.R): each of the nine food groups
# crossed with antibiotics, plus intensity + abx + TPN + EN and the pid / timebin
# random effects; only the food-group columns are replaced by their Z-scores.
# Forest styling matches F2d (shaded abx-interaction bands, blue interaction labels,
# red = 95% CrI clear of zero).

source(here::here("reproduce", "human", "_human_helpers.R"))
if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)
key <- food_key()

# ---- data: F2d diversity data with the food groups Z-scored -----------------
data <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) |>
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid),
         across(starts_with("fg_"), ~ as.numeric(scale(.x))))   # Z-score each food group

# ---- model (priors identical to the F2d diversity family) -------------------
diversity_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
}
fg_vars <- data |> select(starts_with("fg_")) |> colnames()
formula_string <- paste(
  "log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
  paste(paste(fg_vars, "empirical", sep = "*"), collapse = " + "),
  "+ (1 | pid) + (1 | timebin)")

fit <- brm(bf(as.formula(formula_string)), data = data, prior = diversity_priors(),
           warmup = 1000, iter = 3000, chains = 4, cores = 4, seed = 123,
           silent = 2, refresh = 0, control = list(adapt_delta = 0.99),
           backend = brms_backend,
           file = cache_path("E4h_fit_fg_zscored"), file_refit = "on_change")

results_df <- fixef(fit, probs = c(0.025, 0.975)) |> as.data.frame() |>
  rownames_to_column("term") |>
  transmute(effect = "fixed", term, estimate = Estimate,
            conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2))
write_csv(results_df, cache_path("E4h_results_df_fg_zscored.csv"))

# ---- forest plot (F2d styling) ----------------------------------------------
forest_effects <- function(results, label_dict, level_order) {
  results |>
    filter(effect == "fixed") |>
    mutate(clean_term = term |>
        str_replace("empiricalTRUE$", "abx") |>
        str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") |>
        str_replace_all(label_dict) |>
        str_replace("empiricalTRUE:", "abx * ") |>
        str_replace_all("_", " ")) |>
    filter(!str_detect(clean_term, "intensity")) |>
    mutate(is_significant = (conf.low * conf.high) >= 0,
           clean_term = factor(clean_term, levels = level_order))
}
forest_plot <- function(cleaned, band_fill, band_alpha, title) {
  shading <- cleaned |> mutate(y_numeric = as.numeric(clean_term)) |>
    filter(str_detect(clean_term, "\\*"))
  ggplot(cleaned, aes(x = estimate, y = clean_term)) +
    geom_rect(data = shading, aes(ymin = y_numeric - 0.5, ymax = y_numeric + 0.5,
              xmin = -Inf, xmax = Inf), fill = band_fill, alpha = band_alpha,
              inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "solid", color = "blue", linewidth = 0.8) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                    size = 0.25, linewidth = 1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                     str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
    labs(x = "ln (diversity) change", y = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none", axis.text.y = element_markdown(),
          plot.title = element_text(hjust = 0.5, size = 11), aspect.ratio = 1.6)
}

fg_dict <- setNames(key$shortname, key$fg1_name)
fg_order <- rev(c("abx", "EN", "TPN",
  "abx * Sweets", "Sweets", "abx * Grains", "Grains", "abx * Milk", "Milk",
  "abx * Eggs", "Eggs", "abx * Legumes", "Legumes", "abx * Meats", "Meats",
  "abx * Fruits", "Fruits", "abx * Oils", "Oils", "abx * Vegetables", "Vegetables"))

e4h <- forest_plot(forest_effects(results_df, fg_dict, fg_order),
                   band_fill = "#FBEADC", band_alpha = 0.7, title = "FNDDS Z-scored")
save_panel(e4h, "E4h_fndds_zscored_forest.pdf", width = 90, height = 150)

sig <- results_df |> filter((conf.low * conf.high) > 0, str_detect(term, "empirical")) |>
  select(term, estimate, conf.low, conf.high)
message("\nsignificant abx / abx-interaction effects:"); print(sig)
