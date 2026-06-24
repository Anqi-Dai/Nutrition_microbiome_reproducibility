# Fit and cache the F2 diversity models (food-group and macronutrient).
#
# Refactor of the front half of 172_F2_c_to_h.Rmd: the two near-duplicate model
# fits share one prior block and one fitting function, both fits cache to disk via
# brms `file =`, and `sample_prior = TRUE` stores prior draws so the prior/posterior
# predictive diagnostics (12) need no extra refit of the main model. A separate
# prior-only fit (cheap, also cached) backs the prior predictive check.
#
# Outputs (all under intermediate_data/, reused on rerun):
#   172_fit_fg_diversity.rds      food-group model (sample_prior = TRUE)
#   172_fit_macro_diversity.rds   macronutrient model (sample_prior = TRUE)
#   172_fit_fg_prior_only.rds     prior-only fit for the prior predictive check
#   172_results_df_main_fg_diversity.csv / 172_results_df_macro_diversity.csv

source(here::here("reproduce", "human", "_human_helpers.R"))

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

# Intake columns are grams; divide by 100 so the coefficients (and the N(0,1)
# prior) sit on a per-100g scale rather than per-gram.
read_diversity_data <- function() {
  read_csv(released("153_combined_META.csv"), show_col_types = FALSE) %>%
    mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
           pid = factor(pid)) %>%
    mutate(across(c(starts_with("fg_"), starts_with("ave_")), ~ .x / 100))
}

# One prior block, shared by both models. A general N(0,1) on every food/macro
# effect, tighter priors on the clinical covariates, informative priors centring
# each conditioning-intensity intercept near ln(diversity) ~ 2.
diversity_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
}

# Build the diversity formula for a set of predictors (each crossed with abx) and
# fit it, caching the compiled + sampled model. `sample_prior = TRUE` keeps the
# posterior unchanged while also drawing from the prior for diagnostics.
fit_diversity_model <- function(data, predictors, cache_file, sample_prior = TRUE) {
  interaction_terms <- paste(predictors, "empirical", sep = "*")
  formula_string <- paste(
    "log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
    paste(interaction_terms, collapse = " + "),
    "+ (1 | pid) + (1 | timebin)")
  brm(bf(as.formula(formula_string)),
      data = data, prior = diversity_priors(),
      warmup = 1000, iter = 3000, chains = 4, cores = 4,
      seed = 123, silent = 2, refresh = 0,
      control = list(adapt_delta = 0.99),
      sample_prior = sample_prior, backend = brms_backend,
      file = cache_file, file_refit = "on_change")
}

# Coefficient table for the forest panels, with the cohort sizes attached. Built
# from fixef() rather than broom.mixed::tidy(): tidy.brmsfit trips over the
# underscore food-group names plus the extra prior_* draws from sample_prior=TRUE.
# fixef gives the same fixed-effect estimate + 95% CI the panels need.
tidy_results <- function(fit, data) {
  fixef(fit, probs = c(0.025, 0.975)) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    transmute(effect = "fixed", term,
              estimate = Estimate,
              conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2),
              n_samples = nrow(data), n_patients = n_distinct(data$pid))
}

data <- read_diversity_data()
fg_vars <- data %>% select(starts_with("fg_")) %>% colnames()
macro_vars <- c("ave_Fat_g", "ave_Fibers_g", "ave_Sugars_g")

message("Fitting food-group diversity model ...")
fg_fit <- fit_diversity_model(data, fg_vars, cache_path("172_fit_fg_diversity"))

message("Fitting macronutrient diversity model ...")
macro_fit <- fit_diversity_model(data, macro_vars, cache_path("172_fit_macro_diversity"))

message("Fitting prior-only food-group model (for the prior predictive check) ...")
prior_fit <- fit_diversity_model(data, fg_vars, cache_path("172_fit_fg_prior_only"),
                                 sample_prior = "only")

write_csv(tidy_results(fg_fit, data),    cache_path("172_results_df_main_fg_diversity.csv"))
write_csv(tidy_results(macro_fit, data), cache_path("172_results_df_macro_diversity.csv"))
message("Diversity model fits cached to ", intermediate_dir())
