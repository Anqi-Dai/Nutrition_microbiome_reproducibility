# Iterative-subsampling robustness check for the F2 diversity model (E5l).
#
# Refactor of R61. To check that the food-group / antibiotic associations are not
# driven by patients with many stool samples, the model is refit on many random
# subsamples that cap each patient at the cohort median number of samples, and the
# share of runs in which each predictor's 95% CI excludes zero is tallied.
#
# The 50 refits are the expensive part, so the tally caches to intermediate_data/
# and reruns just redraw the bar plot. Panel: E5l.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({ library(future); library(furrr) })

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) %>%
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid)) %>%
  mutate(across(starts_with("fg_"), ~ .x / 100))
key <- food_key()
repl <- setNames(key$shortname, key$fg1_name)

median_samples <- meta %>% count(pid) %>% summarise(median_n = median(n)) %>% pull(median_n)
fg_vars <- meta %>% select(starts_with("fg")) %>% colnames()

# Same formula and priors as the main food-group model.
interactions <- paste(paste(fg_vars, "empirical", sep = "*"), collapse = " + ")
sub_formula <- bf(as.formula(paste(
  "log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
  interactions, "+ (1 | pid) + (1 | timebin)")))
sub_priors <-
  prior(normal(0, 1), class = "b") +
  prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
  prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
  prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
  prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
  prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
  prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")

N_ITERATIONS <- 50
N_CORES <- 8
tally_cache <- cache_path("R61_subsample_tally.rds")

if (file.exists(tally_cache)) {
  final_tally <- readRDS(tally_cache)
} else {
  set.seed(1101)
  plan(multisession, workers = N_CORES)
  cmdstan_dir <- path.expand("~/.cmdstan/cmdstan-2.38.0")
  final_tally <- future_map_dfr(seq_len(N_ITERATIONS), function(i) {
    if (requireNamespace("cmdstanr", quietly = TRUE)) {
      try(cmdstanr::set_cmdstan_path(cmdstan_dir), silent = TRUE)
    }
    subsampled <- meta %>% group_by(pid) %>% slice_sample(n = median_samples) %>% ungroup()
    fit <- brm(sub_formula, data = subsampled, prior = sub_priors,
               iter = 2000, warmup = 1000, chains = 2, cores = 1,
               control = list(adapt_delta = 0.99), silent = 2, refresh = 0,
               backend = "cmdstanr")
    b_vars <- variables(fit)
    fixed_names <- b_vars[startsWith(b_vars, "b_") & !startsWith(b_vars, "b_intensity")]
    fit %>%
      gather_draws(!!!syms(fixed_names)) %>%
      median_qi(.width = 0.95) %>%
      mutate(iteration = i, is_significant = !(.lower < 0 & .upper > 0))
  }, .options = furrr_options(seed = TRUE,
       packages = c("brms", "tidybayes", "dplyr", "stringr", "rlang", "posterior", "cmdstanr")))
  plan(sequential)
  saveRDS(final_tally, tally_cache)
}

# Tally how often each predictor was significant, then clean the names. -------
significance_report <- final_tally %>%
  group_by(.variable) %>%
  summarise(significant_runs = sum(is_significant),
            total_runs = N_ITERATIONS,
            proportion_significant = significant_runs / total_runs, .groups = "drop") %>%
  mutate(clean_variable = .variable %>%
           str_replace("^b_", "") %>%
           str_replace("empiricalTRUE$", "abx") %>%
           str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") %>%
           str_replace_all(repl) %>%
           str_replace("empiricalTRUE:", "abx * ") %>%
           str_replace_all("_", " "))

# E5l: proportion-significant bar plot ----------------------------------------
robustness_plot <- significance_report %>%
  ggplot(aes(x = reorder(clean_variable, proportion_significant), y = proportion_significant)) +
  geom_col(fill = "#0072B2", alpha = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
  labs(title = "Robustness of predictor effects on microbiome diversity",
       subtitle = str_glue("Based on {N_ITERATIONS} random subsamples ",
                           "(capped at {median_samples} samples per patient)"),
       x = "Model predictor",
       y = "Proportion of runs where association was significant (95% CI)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title.position = "plot")

save_panel(robustness_plot, "E5l_robustness_subsampling.pdf", width = 6.0, height = 4.2, units = "in")
message("E5l robustness panel written to results/.")
