# Fit and cache the Extended Fig. E4 diversity models.
#
# Refactor of the model halves of R62_E4.Rmd (E4b) and
# R51_cal_intake_predicting_diversity.Rmd (E4c). These extend the F2 diversity
# family (10_fit_diversity_models.R) to two further predictor sets built from the
# raw diet tracker (DTB) rather than the pre-summarised META food columns:
#   E4b  all five macronutrients (Protein, Fat, Starch, Fibers, Sugars), each
#        crossed with antibiotics, from the 2-day-prior intake window
#   E4c  prior-2-day average caloric intake, fit with and without the
#        abx x calories interaction (the two models compared in the panel)
#
# Starch is derived as Carbohydrates - Fibers - Sugars. The macro window uses the
# join_by window sum/2 of R62; the caloric window uses the two-prior-day mean of
# R51; each is kept faithful to its source so the panels stay line-checkable.
#
# Outputs (under intermediate_data/, reused on rerun):
#   R62_fit_macro5_diversity.rds        five-macronutrient model (E4b)
#   R51_fit_cal_interaction.rds         caloric model with abx interaction (E4c)
#   R51_fit_cal_nointeraction.rds       caloric model without interaction (E4c)
#   R62_results_df_macro5.csv
#   R51_results_df_cal_interaction.csv / R51_results_df_cal_nointeraction.csv

source(here::here("reproduce", "human", "_human_helpers.R"))

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

# Shared prior block (identical to the F2 diversity models): general N(0,1) on the
# food/macro/caloric effects, tighter priors on the clinical covariates, and
# informative intercepts near ln(diversity) ~ 2 per conditioning intensity.
diversity_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
}

# One brms call, shared settings, cached by `file =`.
fit_e4_model <- function(formula_string, data, cache_file) {
  brm(bf(as.formula(formula_string)),
      data = data, prior = diversity_priors(),
      warmup = 1000, iter = 3000, chains = 4, cores = 4,
      seed = 123, silent = 2, refresh = 0,
      control = list(adapt_delta = 0.99),
      backend = brms_backend,
      file = cache_file, file_refit = "on_change")
}

# Coefficient table for the forest panels from fixef() (avoids the broom.mixed
# tidy.brmsfit issue with underscore predictor names); same fixed-effect estimate
# and 95% CI the panels consume.
tidy_results <- function(fit, data) {
  fixef(fit, probs = c(0.025, 0.975)) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    transmute(effect = "fixed", term,
              estimate = Estimate,
              conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2),
              n_samples = nrow(data), n_patients = n_distinct(data$pid))
}

# Diet tracker with starch derived, and the stool-sample frame shared by both
# data preps.
dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE) %>%
  mutate(starch_g = Carbohydrates_g - Fibers_g - Sugars_g)

meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)
stool_samples_df <- meta %>%
  select(pid, sdrt, sampleid, simpson_reciprocal, empirical, intensity, EN, TPN)

# E4b data: five macronutrients on the 2-day-prior window -----------------------
# Daily per-patient macro totals (impute-zero via na.rm sum). rmcorr / the panel
# keep only patients with more than one diet day, so the model shares that subset.
daily_macros <- dtb %>%
  group_by(pid, fdrt) %>%
  summarise(Protein = sum(Protein_g, na.rm = TRUE),
            Fat = sum(Fat_g, na.rm = TRUE),
            Starch = sum(starch_g, na.rm = TRUE),
            Fibers = sum(Fibers_g, na.rm = TRUE),
            Sugars = sum(Sugars_g, na.rm = TRUE),
            .groups = "drop") %>%
  add_count(pid) %>%
  filter(n > 1) %>%
  select(-n)

daily_intake_long <- daily_macros %>%
  pivot_longer(Protein:Sugars, names_to = "macronutrient", values_to = "gram")

macro5_data <- stool_samples_df %>%
  mutate(window_start = sdrt - 2, window_end = sdrt - 1) %>%
  left_join(daily_intake_long,
            by = join_by(pid, window_start <= fdrt, window_end >= fdrt)) %>%
  pivot_wider(id_cols = c(pid, sdrt),
              names_from = macronutrient, values_from = gram,
              values_fn = ~ sum(.x, na.rm = TRUE) / 2,
              values_fill = 0, names_prefix = "avg_intake_") %>%
  mutate(across(starts_with("avg_intake_"), ~ .x / 100)) %>%
  right_join(stool_samples_df, by = c("pid", "sdrt")) %>%
  mutate(timebin = cut_width(sdrt, 7, boundary = 0, closed = "left"),
         intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid))

macro5_vars <- macro5_data %>% select(starts_with("avg_intake_")) %>% colnames()
macro5_formula <- paste(
  "log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
  paste(paste(macro5_vars, "empirical", sep = "*"), collapse = " + "),
  "+ (1 | pid) + (1 | timebin)")

# E4c data: prior-2-day average caloric intake ----------------------------------
# Per-patient daily calories, then the mean of the two days before each stool
# sample (na.rm so one missing prior day still yields an average). Scaled to
# per-1000 kcal so the coefficient reads per 1000 kcal.
daily_calories <- dtb %>%
  group_by(pid, fdrt) %>%
  summarise(total_caloric_intake = sum(Calories_kcal, na.rm = TRUE), .groups = "drop")

cal_data <- stool_samples_df %>%
  mutate(day_prior_1 = sdrt - 1, day_prior_2 = sdrt - 2) %>%
  pivot_longer(c(day_prior_1, day_prior_2), names_to = "day_type", values_to = "fdrt") %>%
  left_join(daily_calories, by = c("pid", "fdrt")) %>%
  group_by(pid, sdrt) %>%
  summarize(avg_cal_intake_2_day = mean(total_caloric_intake, na.rm = TRUE), .groups = "drop") %>%
  right_join(stool_samples_df, by = c("pid", "sdrt")) %>%
  mutate(timebin = cut_width(sdrt, 7, boundary = 0, closed = "left"),
         intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid)) %>%
  mutate(across(starts_with("avg_cal_intake_2_day"), ~ .x / 1000))

cal_formula_interaction <- paste(
  "log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
  "avg_cal_intake_2_day*empirical + (1 | pid) + (1 | timebin)")
cal_formula_nointeraction <- paste(
  "log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
  "avg_cal_intake_2_day + empirical + (1 | pid) + (1 | timebin)")

# Fit + cache -------------------------------------------------------------------
message("Fitting five-macronutrient diversity model (E4b) ...")
macro5_fit <- fit_e4_model(macro5_formula, macro5_data, cache_path("R62_fit_macro5_diversity"))

message("Fitting caloric-intake diversity model with interaction (E4c) ...")
cal_fit_int <- fit_e4_model(cal_formula_interaction, cal_data, cache_path("R51_fit_cal_interaction"))

message("Fitting caloric-intake diversity model without interaction (E4c) ...")
cal_fit_noint <- fit_e4_model(cal_formula_nointeraction, cal_data, cache_path("R51_fit_cal_nointeraction"))

write_csv(tidy_results(macro5_fit, macro5_data),  cache_path("R62_results_df_macro5.csv"))
write_csv(tidy_results(cal_fit_int, cal_data),    cache_path("R51_results_df_cal_interaction.csv"))
write_csv(tidy_results(cal_fit_noint, cal_data),  cache_path("R51_results_df_cal_nointeraction.csv"))

# E4j data + fit: does the abx x sweets effect depend on conditioning intensity?
# Refactor of R36: a three-way fg_sweets x abx x intensity interaction model.
# Faithful to R36 and distinct from the F2/E4 fits above, so it does not reuse
# diversity_priors(): intensity is left at its default (alphabetical) factor levels
# so the ablative level is the reference, intensity enters only through the
# interaction expansion (no clean intercept), and the prior block therefore omits
# the intensity-intercept priors. The other food groups stay as additive controls.
e4j_data <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) %>%
  mutate(pid = factor(pid)) %>%
  mutate(across(starts_with("fg_"), ~ .x / 100))

e4j_other_fg <- e4j_data %>% select(starts_with("fg")) %>% select(-fg_sweets) %>% colnames()
e4j_formula <- paste(
  "log(simpson_reciprocal) ~ 0 + TPN + EN +",
  paste(e4j_other_fg, collapse = " + "),
  "+ fg_sweets * empirical * intensity + (1 | pid) + (1 | timebin)")

e4j_priors <- prior(normal(0, 1), class = "b") +
  prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
  prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
  prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE")

message("Fitting three-way sweets x abx x intensity interaction model (E4j) ...")
e4j_fit <- brm(bf(as.formula(e4j_formula)),
               data = e4j_data, prior = e4j_priors,
               warmup = 1000, iter = 3000, chains = 4, cores = 4,
               seed = 123, silent = 2, refresh = 0,
               backend = brms_backend,
               file = cache_path("R36_fit_sweets_intensity"), file_refit = "on_change")

# The panel is the marginal sweets slope within antibiotic-exposed patients for
# each intensity: the baseline interaction is the (reference) ablative group, and
# each non-reference level adds its three-way term. Cache the conditional draws so
# 17 stays plot-only.
e4j_draws <- as_draws_df(e4j_fit) %>%
  transmute(
    Myeloablative = `b_fg_sweets:empiricalTRUE`,
    `Reduced Intensity` = `b_fg_sweets:empiricalTRUE` + `b_fg_sweets:empiricalTRUE:intensityreduced`,
    Nonmyeloablative = `b_fg_sweets:empiricalTRUE` + `b_fg_sweets:empiricalTRUE:intensitynonablative`)
write_csv(e4j_draws, cache_path("R36_e4j_conditional_draws.csv"))

# E4e data + fit: split the FNDDS group-9 "sweets" lump into its two-digit subgroups
# Refactor of R27_sub_groups_of_sweets.Rmd. The F2/E4b food-group model collapses
# all of FNDDS major group 9 (Sugars, Sweets, and Beverages) into one "sweets"
# predictor; here group 9 is split on the two-digit FNDDS code into its subgroups
# (91 sugars & sweets, 92 non-alcoholic beverages, 94 water, 95 formulated
# beverages; 93 alcoholic beverages never appears in the cohort), each kept as its
# own predictor crossed with antibiotics, while groups 1-8 stay as the usual eight
# food groups. Same 2-day-prior intake window and shared diversity prior block as
# E4b, so the abx interactions stay comparable across panels.
fndds_daily_intake <- dtb %>%
  mutate(Food_code = as.character(Food_code)) %>%
  mutate(food_category = case_when(
    str_starts(Food_code, "91") ~ "Sugars_and_sweets",
    str_starts(Food_code, "92") ~ "Nonalcoholic_beverages",
    str_starts(Food_code, "93") ~ "Alcoholic_beverages",
    str_starts(Food_code, "94") ~ "Water",
    str_starts(Food_code, "95") ~ "Formulated_beverages",
    str_starts(Food_code, "1") ~ "fg_milk",
    str_starts(Food_code, "2") ~ "fg_meat",
    str_starts(Food_code, "3") ~ "fg_egg",
    str_starts(Food_code, "4") ~ "fg_legume",
    str_starts(Food_code, "5") ~ "fg_grain",
    str_starts(Food_code, "6") ~ "fg_fruit",
    str_starts(Food_code, "7") ~ "fg_veggie",
    str_starts(Food_code, "8") ~ "fg_oils",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(food_category)) %>%
  group_by(pid, fdrt, food_category) %>%
  summarise(total_dehydrated_weight = sum(dehydrated_weight, na.rm = TRUE),
            .groups = "drop")

fndds_data <- stool_samples_df %>%
  mutate(window_start = sdrt - 2, window_end = sdrt - 1) %>%
  left_join(fndds_daily_intake,
            by = join_by(pid, window_start <= fdrt, window_end >= fdrt)) %>%
  pivot_wider(id_cols = c(pid, sdrt),
              names_from = food_category, values_from = total_dehydrated_weight,
              values_fn = ~ sum(.x, na.rm = TRUE) / 2,
              values_fill = 0, names_prefix = "avg_intake_") %>%
  mutate(across(starts_with("avg_intake_"), ~ .x / 100)) %>%
  right_join(stool_samples_df, by = c("pid", "sdrt")) %>%
  mutate(timebin = cut_width(sdrt, 7, boundary = 0, closed = "left"),
         intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid))

fndds_vars <- fndds_data %>% select(starts_with("avg_intake_")) %>% colnames()
fndds_formula <- paste(
  "log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
  paste(paste(fndds_vars, "empirical", sep = "*"), collapse = " + "),
  "+ (1 | pid) + (1 | timebin)")

message("Fitting FNDDS two-digit subgroup diversity model (E4e) ...")
fndds_fit <- fit_e4_model(fndds_formula, fndds_data, cache_path("R27_fit_fndds_subgroups"))
write_csv(tidy_results(fndds_fit, fndds_data), cache_path("R27_results_df_fndds.csv"))

# Per-subgroup exposure counts for the nomenclature tree panel (E4f): how many
# stool samples (and distinct patients) had any intake of each FNDDS subgroup in
# the 2-day-prior window. Deterministic, but written here so the tree panel reads
# the same model frame the forest is fit on (93 alcoholic beverages is absent from
# the cohort, so it reports 0 / 0).
fndds_exposure <- fndds_data %>%
  select(pid, starts_with("avg_intake_")) %>%
  pivot_longer(-pid, names_to = "food_category", values_to = "avg_intake",
               names_prefix = "avg_intake_") %>%
  group_by(food_category) %>%
  summarise(samples = sum(avg_intake > 0),
            patients = n_distinct(pid[avg_intake > 0]), .groups = "drop") %>%
  filter(food_category %in% c("Sugars_and_sweets", "Nonalcoholic_beverages",
                              "Water", "Formulated_beverages"))
write_csv(fndds_exposure, cache_path("R27_fndds_exposure_counts.csv"))

message("E4 model fits cached to ", intermediate_dir())
