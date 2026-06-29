# A4 three-panel forest figure on the prior-5-day diet exposure window.
#
# Cohort: a stool sample qualifies if EACH of its 5 prior days is covered -- either
# at least one dietary entry, or a documented zero-intake day (072). That is the
# 801-sample / 146-patient set. Diet exposure for each food group / macronutrient
# is the prior-5-day average (sum over [sdrt-5, sdrt-1] / 5; zero-eating days carry
# 0). Models, priors and structure are identical to the published F2/F4b fits.
#
# Panels (forest plots, abx-interaction rows shaded, red = CI clear of zero):
#   A  F2d-style food-group diversity effects (log Simpson)
#   B  F4b-style E. faecium (asv_1) CLR effects
#   C  F2e-style macronutrient diversity effects (log Simpson)
#
# Output: results/F2_F4b_forest_5day_A4.pdf

source(here::here("reproduce", "human", "_human_helpers.R"))
if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

WINDOW <- 5L
key <- food_key()

# ---- priors / fit helpers (identical to 10/16 and 40) ---------------------
diversity_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
}
taxon_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE")
}
fit_brm <- function(formula_string, data, priors, cache_file) {
  brm(bf(as.formula(formula_string)), data = data, prior = priors,
      warmup = 1000, iter = 3000, chains = 4, cores = 4, seed = 123,
      silent = 2, refresh = 0, control = list(adapt_delta = 0.99),
      backend = brms_backend, file = cache_file, file_refit = "on_change")
}
tidy_results <- function(fit, data) {
  fixef(fit, probs = c(0.025, 0.975)) |> as.data.frame() |>
    rownames_to_column("term") |>
    transmute(effect = "fixed", term, estimate = Estimate,
              conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2))
}

# ---- cohort: all 5 prior days covered (entry or documented zero) -----------
dtb  <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE)
meta <- read_csv(released("R59_meta_expanded.csv"), show_col_types = FALSE)
zero <- read_csv(released("072_total_patients_zero_eating_days_pid.csv"),
                 show_col_types = FALSE) |> distinct(pid, fdrt)

valid_days <- bind_rows(dtb |> distinct(pid, fdrt), zero) |>
  distinct(pid, fdrt) |> mutate(has = TRUE)

cohort <- meta |> select(pid, sdrt, sampleid) |>
  tidyr::crossing(k = 1:WINDOW) |> mutate(fdrt = sdrt - k) |>
  left_join(valid_days, by = c("pid", "fdrt")) |>
  group_by(sampleid) |> summarise(nok = sum(!is.na(has)), .groups = "drop") |>
  filter(nok == WINDOW) |> pull(sampleid)

samples <- meta |> filter(sampleid %in% cohort) |> select(pid, sdrt, sampleid)
message("Cohort (all ", WINDOW, " prior days covered): ", length(cohort),
        " samples, ", n_distinct(samples$pid), " patients")

# ---- prior-5-day diet exposure (sum / WINDOW) -----------------------------
fg_map <- c("1" = "fg_milk", "2" = "fg_meat", "3" = "fg_egg", "4" = "fg_legume",
            "5" = "fg_grain", "6" = "fg_fruit", "7" = "fg_veggie",
            "8" = "fg_oils", "9" = "fg_sweets")
daily_fg <- dtb |>
  mutate(fg = fg_map[substr(as.character(Food_code), 1, 1)]) |>
  filter(!is.na(fg)) |>
  group_by(pid, fdrt, fg) |> summarise(w = sum(dehydrated_weight), .groups = "drop")
daily_macro <- dtb |> group_by(pid, fdrt) |>
  summarise(across(c(Fat_g, Fibers_g, Sugars_g), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

fg5 <- samples |> mutate(ws = sdrt - WINDOW, we = sdrt - 1) |>
  left_join(daily_fg, by = join_by(pid, ws <= fdrt, we >= fdrt)) |>
  filter(!is.na(fg)) |>
  pivot_wider(id_cols = sampleid, names_from = fg, values_from = w,
              values_fn = ~ sum(.x, na.rm = TRUE) / WINDOW, values_fill = 0)
macro5 <- samples |> mutate(ws = sdrt - WINDOW, we = sdrt - 1) |>
  left_join(daily_macro, by = join_by(pid, ws <= fdrt, we >= fdrt)) |>
  group_by(sampleid) |>
  summarise(ave_Fat_g = sum(Fat_g, na.rm = TRUE) / WINDOW,
            ave_Fibers_g = sum(Fibers_g, na.rm = TRUE) / WINDOW,
            ave_Sugars_g = sum(Sugars_g, na.rm = TRUE) / WINDOW, .groups = "drop")

model_data <- meta |> filter(sampleid %in% cohort) |>
  select(pid, sdrt, sampleid, simpson_reciprocal, empirical, intensity,
         EN, TPN, timebin, asv_1_clr) |>
  left_join(fg5, by = "sampleid") |> left_join(macro5, by = "sampleid") |>
  mutate(across(c(starts_with("fg_"), starts_with("ave_")), ~ replace_na(.x, 0) / 100),
         intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid))

n_s <- nrow(model_data); n_p <- n_distinct(model_data$pid)

# ---- fit the three models -------------------------------------------------
fg_vars    <- grep("^fg_", names(model_data), value = TRUE)
macro_vars <- c("ave_Fat_g", "ave_Fibers_g", "ave_Sugars_g")
rhs <- function(v) paste("0 + intensity + empirical + TPN + EN +",
  paste(paste(v, "empirical", sep = "*"), collapse = " + "), "+ (1 | pid) + (1 | timebin)")

message("Fitting food-group diversity ...")
fg_fit <- fit_brm(paste("log(simpson_reciprocal) ~", rhs(fg_vars)), model_data,
                  diversity_priors(), cache_path("172_fit_fg_diversity_5day801"))
message("Fitting macronutrient diversity ...")
macro_fit <- fit_brm(paste("log(simpson_reciprocal) ~", rhs(macro_vars)), model_data,
                     diversity_priors(), cache_path("172_fit_macro_diversity_5day801"))
message("Fitting E. faecium asv_1 CLR ...")
asv1_fit <- fit_brm(paste("asv_1_clr ~", rhs(fg_vars)), model_data,
                    taxon_priors(), cache_path("R63_fit_asv1_5day801"))

# ---- forest builders (match 11 / 41 styling) ------------------------------
forest_clean <- function(results, label_dict, level_order) {
  results |> filter(effect == "fixed") |>
    mutate(clean_term = term |>
             str_replace("empiricalTRUE$", "abx") |>
             str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") |>
             str_replace_all(label_dict) |>
             str_replace("empiricalTRUE:", "abx * ") |>
             str_replace_all("_", " ")) |>
    filter(!str_detect(clean_term, "intensity")) |>
    mutate(is_significant = (conf.low * conf.high) > 0,
           clean_term = factor(clean_term, levels = level_order))
}
forest_plot <- function(cleaned, xlab, title, band_fill, band_alpha, title_italic = FALSE) {
  shading <- cleaned |> mutate(y_numeric = as.numeric(clean_term)) |>
    filter(str_detect(clean_term, "\\*"))
  ggplot(cleaned, aes(x = estimate, y = clean_term)) +
    geom_rect(data = shading, aes(ymin = y_numeric - 0.5, ymax = y_numeric + 0.5,
              xmin = -Inf, xmax = Inf), fill = band_fill, alpha = band_alpha,
              inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "solid", color = "blue", linewidth = 0.8) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                    size = 0.25, linewidth = 0.9) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                     str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
    labs(x = xlab, y = "", title = title) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none", axis.text.y = element_markdown(),
          plot.title = element_text(hjust = 0.5, size = 11,
                                    face = if (title_italic) "italic" else "plain"),
          axis.text = element_text(size = 9), axis.title = element_text(size = 10))
}

fg_dict <- setNames(key$shortname, key$fg1_name)
fg_order <- rev(c("abx", "EN", "TPN",
  "abx * Sweets", "Sweets", "abx * Grains", "Grains", "abx * Milk", "Milk",
  "abx * Eggs", "Eggs", "abx * Legumes", "Legumes", "abx * Meats", "Meats",
  "abx * Fruits", "Fruits", "abx * Oils", "Oils", "abx * Vegetables", "Vegetables"))
macro_dict <- c("ave_Fat_g" = "Fat", "ave_Fibers_g" = "Fiber", "ave_Sugars_g" = "Sugars")
macro_order <- rev(c("abx", "EN", "TPN",
  "abx * Sugars", "Sugars", "abx * Fiber", "Fiber", "abx * Fat", "Fat"))

pA <- forest_plot(forest_clean(tidy_results(fg_fit, model_data), fg_dict, fg_order),
                  "ln(diversity) change", "Food groups -> microbiome diversity",
                  "#FBEADC", 0.7)
pB <- forest_plot(forest_clean(tidy_results(asv1_fit, model_data), fg_dict, fg_order),
                  "CLR E. faecium change", "Food groups -> E. faecium",
                  "#FBEADC", 0.7, title_italic = FALSE)
pC <- forest_plot(forest_clean(tidy_results(macro_fit, model_data), macro_dict, macro_order),
                  "ln(diversity) change", "Macronutrients -> microbiome diversity",
                  "#6B8E2340", 0.3)

# ---- assemble A4 (210 x 297 mm) -------------------------------------------
# Layout: a (food-group diversity) top-left, c (macronutrient diversity) below it
# in the SAME left column with their x-axes aligned (align = "v", axis = "lr");
# b (E. faecium) top-right. Panel heights are proportional to row counts (a/b = 21
# rows, c = 9) so every coefficient row is the same height across all three panels.
# Page margins via the draw_plot inset.
pad <- function(p) p + theme(plot.margin = margin(3, 6, 3, 6, "mm"))
n_fg <- 21; n_mac <- 9   # coefficient rows per panel (for equal row height)

left_col  <- cowplot::plot_grid(pad(pA), pad(pC), ncol = 1,
                                rel_heights = c(n_fg, n_mac), labels = c("a", "c"),
                                align = "v", axis = "lr")
right_col <- cowplot::plot_grid(pad(pB), NULL, ncol = 1,
                                rel_heights = c(n_fg, n_mac), labels = c("b", ""))
body <- cowplot::plot_grid(left_col, right_col, ncol = 2)

title <- cowplot::ggdraw() +
  cowplot::draw_label(
    sprintf(paste0("Prior-5-day diet exposure: n = %d stool samples, %d patients\n",
                   "(each of the 5 prior days had >=1 dietary entry or a confirmed zero-intake day)"),
            n_s, n_p), size = 11, x = 0.5, hjust = 0.5, lineheight = 1.1)
withtitle <- cowplot::plot_grid(title, body, ncol = 1, rel_heights = c(0.05, 1))

# page margins: inset the whole figure inside the A4 sheet
fig <- cowplot::ggdraw() +
  cowplot::draw_plot(withtitle, x = 0.045, y = 0.035, width = 0.91, height = 0.93)

out_dir <- here::here("results"); if (!dir.exists(out_dir)) dir.create(out_dir)
ggsave(file.path(out_dir, "F2_F4b_forest_5day_A4.pdf"), fig,
       width = 210, height = 297, units = "mm")
message("wrote results/F2_F4b_forest_5day_A4.pdf  (n=", n_s, " samples, ", n_p, " patients)")
