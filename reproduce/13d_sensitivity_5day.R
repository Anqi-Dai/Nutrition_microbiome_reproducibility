# Sensitivity analysis: do the F2 diversity (food-group, macronutrient) and F4b
# E. faecium (asv_1) results hold if diet exposure is summarised over the prior
# FIVE days instead of the published two?
#
# The published models use the 2-day-prior diet exposure carried in the released
# meta (fg_* / ave_*, each = sum over [sdrt-2, sdrt-1] / 2). Here we recompute the
# same exposures over a 5-day window (sum over [sdrt-5, sdrt-1] / 5) directly from
# 152_combined_DTB, keeping everything else identical (same 1009-sample cohort,
# same covariates, priors, and model structure). Documented zero-eating days carry
# zero intake; with the fixed /W denominator they contribute 0 just like the
# 2-day window, so "counting zero-eating days" needs no special handling here --
# the window simply spans more days.
#
# Outputs (intermediate_data/, *_5day):
#   172_fit_fg_diversity_5day.rds / _results_*    F2 food-group, 5-day exposure
#   172_fit_macro_diversity_5day.rds / _results_* F2 macronutrient, 5-day exposure
#   R63_fit_asv1_5day.rds / _results_*            F4b E. faecium CLR, 5-day exposure
# and a side-by-side comparison of the key effects vs the published 2-day fits.

source(here::here("reproduce", "human", "_human_helpers.R"))
if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

WINDOW <- as.integer(Sys.getenv("DIET_WINDOW", "5"))   # prior-day window length
message("Diet exposure window: prior ", WINDOW, " days")

# ---- priors (identical to 10/16 for diversity, 40 for taxon) ---------------
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
    transmute(term, estimate = round(Estimate, 2),
              conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2),
              sig = if_else(conf.low > 0 | conf.high < 0, "*", ""),
              n_samples = nrow(data), n_patients = n_distinct(data$pid))
}

# ---- recompute diet exposure over the W-day prior window from 152 ----------
dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE)
meta <- read_csv(released("R59_meta_expanded.csv"), show_col_types = FALSE)
samples <- meta |> select(pid, sdrt, sampleid)

fg_map <- c("1" = "fg_milk", "2" = "fg_meat", "3" = "fg_egg", "4" = "fg_legume",
            "5" = "fg_grain", "6" = "fg_fruit", "7" = "fg_veggie",
            "8" = "fg_oils", "9" = "fg_sweets")

# daily per-patient food-group dehydrated weight (digit-1 of Food_code -> group)
daily_fg <- dtb |>
  mutate(fg = fg_map[substr(as.character(Food_code), 1, 1)]) |>
  filter(!is.na(fg)) |>
  group_by(pid, fdrt, fg) |>
  summarise(w = sum(dehydrated_weight), .groups = "drop")

# daily per-patient macronutrient totals
daily_macro <- dtb |>
  group_by(pid, fdrt) |>
  summarise(across(c(Fat_g, Fibers_g, Sugars_g), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

# window sum / W; samples with no diet in the window get 0 (zero-eating or absent)
window_wide <- function(daily, value_cols, name_col = NULL) {
  j <- samples |> mutate(ws = sdrt - WINDOW, we = sdrt - 1) |>
    left_join(daily, by = join_by(pid, ws <= fdrt, we >= fdrt))
  if (!is.null(name_col)) {
    j <- j |> filter(!is.na(.data[[name_col]])) |>
      pivot_wider(id_cols = sampleid, names_from = all_of(name_col),
                  values_from = all_of(value_cols),
                  values_fn = ~ sum(.x, na.rm = TRUE) / WINDOW, values_fill = 0)
  } else {
    j <- j |> group_by(sampleid) |>
      summarise(across(all_of(value_cols), ~ sum(.x, na.rm = TRUE) / WINDOW),
                .groups = "drop")
  }
  tibble(sampleid = samples$sampleid) |> left_join(j, by = "sampleid") |>
    mutate(across(-sampleid, ~ replace_na(.x, 0)))
}

fg5    <- window_wide(daily_fg, "w", name_col = "fg")
macro5 <- window_wide(daily_macro, c("Fat_g", "Fibers_g", "Sugars_g")) |>
  rename(ave_Fat_g = Fat_g, ave_Fibers_g = Fibers_g, ave_Sugars_g = Sugars_g)

# ---- assemble model data: covariates from meta, exposures from the W-day window
base <- meta |>
  select(pid, sdrt, sampleid, simpson_reciprocal, empirical, intensity,
         EN, TPN, timebin, asv_1_clr) |>
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid))

model_data <- base |>
  left_join(fg5, by = "sampleid") |>
  left_join(macro5, by = "sampleid") |>
  mutate(across(c(starts_with("fg_"), starts_with("ave_")), ~ .x / 100))

fg_vars    <- grep("^fg_", names(model_data), value = TRUE)
macro_vars <- c("ave_Fat_g", "ave_Fibers_g", "ave_Sugars_g")
rhs <- function(vars) paste(
  "0 + intensity + empirical + TPN + EN +",
  paste(paste(vars, "empirical", sep = "*"), collapse = " + "),
  "+ (1 | pid) + (1 | timebin)")

# ---- fit the three models with 5-day exposure -----------------------------
message("Fitting F2 food-group diversity (", WINDOW, "-day) ...")
fg_fit <- fit_brm(paste("log(simpson_reciprocal) ~", rhs(fg_vars)),
                  model_data, diversity_priors(),
                  cache_path(paste0("172_fit_fg_diversity_", WINDOW, "day")))

message("Fitting F2 macronutrient diversity (", WINDOW, "-day) ...")
macro_fit <- fit_brm(paste("log(simpson_reciprocal) ~", rhs(macro_vars)),
                     model_data, diversity_priors(),
                     cache_path(paste0("172_fit_macro_diversity_", WINDOW, "day")))

message("Fitting F4b E. faecium asv_1 CLR (", WINDOW, "-day) ...")
asv1_fit <- fit_brm(paste("asv_1_clr ~", rhs(fg_vars)),
                    model_data, taxon_priors(),
                    cache_path(paste0("R63_fit_asv1_", WINDOW, "day")))

fg_res    <- tidy_results(fg_fit, model_data)
macro_res <- tidy_results(macro_fit, model_data)
asv1_res  <- tidy_results(asv1_fit, model_data)
write_csv(fg_res,    cache_path(paste0("172_results_fg_diversity_", WINDOW, "day.csv")))
write_csv(macro_res, cache_path(paste0("172_results_macro_diversity_", WINDOW, "day.csv")))
write_csv(asv1_res,  cache_path(paste0("R63_results_asv1_", WINDOW, "day.csv")))

# ---- side-by-side vs the published 2-day fits ------------------------------
pub_path <- function(f) cache_path(f)
compare <- function(new_res, pub_file, label, terms) {
  pub <- read_csv(pub_path(pub_file), show_col_types = FALSE) |>
    transmute(term, est_2day = estimate, lo_2day = conf.low, hi_2day = conf.high)
  out <- new_res |>
    transmute(term, est_5day = estimate, lo_5day = conf.low, hi_5day = conf.high, sig_5day = sig) |>
    left_join(pub, by = "term") |>
    filter(term %in% terms)
  cat("\n==== ", label, " (2-day published vs 5-day) ====\n"); print(as.data.frame(out))
}

compare(fg_res, "172_results_df_main_fg_diversity.csv", "F2 food-group diversity",
        c("empiricalTRUE", "empiricalTRUE:fg_sweets", "empiricalTRUE:fg_grain",
          "fg_sweets", "fg_grain"))
compare(macro_res, "172_results_df_macro_diversity.csv", "F2 macronutrient diversity",
        c("empiricalTRUE", "empiricalTRUE:ave_Sugars_g", "empiricalTRUE:ave_Fibers_g",
          "ave_Sugars_g", "ave_Fibers_g"))
compare(asv1_res, "R63_results_df_asv1.csv", "F4b E. faecium asv_1 CLR",
        c("empiricalTRUE", "empiricalTRUE:fg_sweets", "fg_sweets"))

message("\nDone. Full 5-day coefficient tables written to intermediate_data/*_", WINDOW, "day.csv")
