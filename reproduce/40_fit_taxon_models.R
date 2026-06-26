# Fit and cache the human taxon-abundance models (Enterococcus / genus CLR family).
#
# Refactor of R63_Enterococcus_asv1_outcome.Rmd. Two model families share one
# right-hand side (the F2 diversity predictor set: conditioning intensity,
# antibiotics, TPN/EN and the nine food groups each crossed with antibiotics):
#   F4b   one model with E. faecium (asv_1) CLR as outcome
#   E7a   one model per prevalent genus (CLR as outcome), summarised into the
#         heatmap of food-group effects across genera
#
# These differ from the diversity family in the outcome (a centred-log-ratio
# abundance, not log Simpson) and in the prior block: intensity enters via
# `0 + intensity` but keeps the general N(0,1) prior (no informative diversity
# intercepts), so this does not reuse 16's diversity_priors().
#
# Outputs (under intermediate_data/, reused on rerun):
#   R63_fit_asv1.rds                     E. faecium CLR model (F4b, E7e source)
#   R63_results_df_asv1.csv              its fixed-effect table (F4b forest)
#   R63_genus_clr_all_models_results.csv per-genus fixed-effect summaries (E7a)

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages(library(furrr))

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

# Shared prior block. General N(0,1) on every coefficient (intensity intercepts
# included), tighter priors on the clinical covariates and a moderately
# informative prior on the antibiotics main effect.
taxon_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE")
}

# One brms call, shared sampler settings.
fit_taxon_brm <- function(formula_string, data, cache_file = NULL, cores = 4) {
  brm(bf(as.formula(formula_string)),
      data = data, prior = taxon_priors(),
      warmup = 1000, iter = 3000, chains = 4, cores = cores,
      seed = 123, silent = 2, refresh = 0,
      control = list(adapt_delta = 0.99),
      backend = brms_backend,
      file = cache_file, file_refit = "on_change")
}

# Meta with the model predictors. fg_* ship in grams; the models use per-100g
# units (the slope then reads per 100 g eaten). asv_1_clr is the E. faecium CLR
# carried in the released meta.
meta <- read_csv(released("R59_meta_expanded.csv"), show_col_types = FALSE) |>
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid),
         across(starts_with("fg_"), ~ .x / 100))

all_food_vars <- meta |> select(starts_with("fg")) |> colnames()
interaction_terms <- paste(all_food_vars, "empirical", sep = "*")
base_rhs <- paste(
  "0 + intensity + empirical + TPN + EN +",
  paste(interaction_terms, collapse = " + "),
  "+ (1 | pid) + (1 | timebin)")

# E. faecium (asv_1) model (F4b) ------------------------------------------------
message("Fitting E. faecium (asv_1) CLR model (F4b) ...")
asv1_formula <- paste("asv_1_clr ~", base_rhs)
asv1_fit <- fit_taxon_brm(asv1_formula, meta, cache_path("R63_fit_asv1"))

# Fixed-effect table for the forest, from fixef() (avoids the broom.mixed
# tidy.brmsfit issue with underscore predictor names). Same estimate + 95% CI
# the panel consumes.
asv1_results <- fixef(asv1_fit, probs = c(0.025, 0.975)) |>
  as.data.frame() |>
  rownames_to_column("term") |>
  transmute(effect = "fixed", term,
            estimate = Estimate,
            conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2),
            n_samples = nrow(meta), n_patients = n_distinct(meta$pid))
write_csv(asv1_results, cache_path("R63_results_df_asv1.csv"))

# Effect size of the abx * sweets term, on the CLR and fold scale (panel text).
asv1_results |>
  filter(term == "empiricalTRUE:fg_sweets") |>
  mutate(across(c(estimate, conf.low, conf.high), list(fold = ~ round(exp(.x), 2)),
                .names = "{.col}_fold")) |>
  print()

# Per-genus models (E7a) --------------------------------------------------------
# Keep the prevalent, reasonably abundant genera likely to exist in the human
# gut: present (relab > 0.002) in more than 10% of samples, dropping two
# environmental genera, then join each genus CLR onto the meta as a candidate
# outcome.
asv_relab_genus <- read_csv(released("171_quality_asv_relab_pident97_genus.csv"), show_col_types = FALSE)
clr_res <- read_csv(released("171_genus_CLR_res.csv"), show_col_types = FALSE)

df_relab <- asv_relab_genus |>
  group_by(sampleid, genus) |>
  summarize(relab = sum(count_relative, na.rm = TRUE), .groups = "drop") |>
  filter(sampleid %in% meta$sampleid)

target_genera <- df_relab |>
  group_by(genus) |>
  count(relab > 0.002) |>
  rename(criteria = 2) |>
  filter(criteria == "TRUE") |>
  arrange(-n) |>
  filter(!is.na(genus)) |>
  filter(!genus %in% c("Ruthenibacterium", "Drancourtella")) |>
  mutate(perc = round(n / nrow(meta) * 100, 0)) |>
  filter(perc > 10) |>
  pull(genus)

clr_wide <- clr_res |>
  filter(genus %in% target_genera) |>
  spread("genus", "clr")

meta_genus_clr <- meta |> inner_join(clr_wide, by = "sampleid")

taxa_outcomes <- meta_genus_clr |> select(all_of(sort(target_genera))) |> colnames()
message("Fitting ", length(taxa_outcomes), " per-genus CLR models (E7a) ...")

# Per-genus fit. Returns the fixed-effect summary at tiered credible intervals
# (94/97/99%) plus the fraction of posterior draws strictly positive (used for
# the dendrogram ordering). Workers run as fresh sessions, so re-resolve the
# CmdStan path inside the function.
fit_genus_model <- function(outcome, data, rhs, backend) {
  suppressPackageStartupMessages(library(brms))
  if (backend == "cmdstanr") {
    if (is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL))) {
      try(cmdstanr::set_cmdstan_path(path.expand("~/.cmdstan/cmdstan-2.38.0")), silent = TRUE)
    }
  }
  model_formula <- brms::bf(as.formula(paste(outcome, "~", rhs)))
  priors <- brms::prior(normal(0, 1), class = "b") +
    brms::prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    brms::prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    brms::prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE")

  fit <- brms::brm(formula = model_formula, data = data, prior = priors,
                   warmup = 1000, iter = 3000, chains = 4, cores = 1,
                   seed = 123, silent = 2, refresh = 0,
                   control = list(adapt_delta = 0.99), backend = backend)

  ci_probs <- c(0.005, 0.015, 0.025, 0.03, 0.97, 0.975, 0.985, 0.995)
  results_summary <- brms::fixef(fit, probs = ci_probs) |>
    tibble::as_tibble(rownames = "term") |>
    dplyr::mutate(outcome = outcome, .before = 1)

  prob_positive <- brms::as_draws_df(fit) |>
    dplyr::select(dplyr::starts_with("b_")) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), \(x) mean(x > 0))) |>
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "term",
                        values_to = "fraction_positive") |>
    dplyr::mutate(term = stringr::str_remove(term, "^b_"))

  dplyr::left_join(results_summary, prob_positive, by = "term")
}

n_workers <- as.integer(Sys.getenv("GENUS_WORKERS", unset = "6"))
if (n_workers > 1) {
  future::plan(multisession, workers = n_workers)
  genus_results <- taxa_outcomes |>
    set_names() |>
    future_map_dfr(\(o) fit_genus_model(o, meta_genus_clr, base_rhs, brms_backend),
                   .options = furrr_options(seed = TRUE), .progress = TRUE)
  future::plan(sequential)
} else {
  genus_results <- taxa_outcomes |>
    set_names() |>
    map_dfr(\(o) fit_genus_model(o, meta_genus_clr, base_rhs, brms_backend))
}

write_csv(genus_results, cache_path("R63_genus_clr_all_models_results.csv"))
message("Taxon model fits cached to ", intermediate_dir())
