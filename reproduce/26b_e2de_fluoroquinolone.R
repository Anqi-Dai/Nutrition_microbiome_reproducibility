# Extended Fig. E2 d,e: the simplified diversity model that carries an explicit
# term for prophylactic fluoroquinolone exposure.
#
# There is no reference .Rmd for this pair of panels (the fit was never scripted in
# this repo), so the model is rebuilt from the formula printed inside the published
# panel d and from the E2 legend:
#
#   ln(alpha-diversity) ~ conditioning + total parenteral nutrition +
#                         enteral nutrition + broad-spectrum antibiotics +
#                         prophylactic fluoroquinolones + (1|patient_id) + (1|week)
#
# It is the F2d food-group model with every dietary term (and therefore every
# antibiotic x diet interaction) dropped, and one new binary predictor added. So it
# reuses the F2 machinery verbatim: the same cohort (153_combined_META.csv, 1009
# samples / 158 patients), the same conditioning-intensity cell-mean intercepts
# (0 + intensity), the same varying intercepts for patient and for week relative to
# transplant (timebin), and the same shared prior block, diversity_priors(). The new
# fluoroquinolone coefficient falls under that block's general N(0,1) on class "b".
#
#   E2d  forest of the four exposure coefficients (posterior mean + 95% CrI),
#        red where the interval excludes zero
#   E2e  the varying-intercept coefficients of the week term, one per week bin
#
# Nothing about broad-spectrum exposure changes: `empirical` in the released meta
# table is exactly "any broad-spectrum drug in the two days before the sample", and
# the script asserts that against Data_S4 before fitting.
#
# Outputs: intermediate_data/E2de_fit_fluoroquinolone.rds (cached fit),
#          intermediate_data/E2de_results_df_fluoroquinolone.csv,
#          results/E2d_fluoroquinolone_forest.pdf, results/E2e_week_intercepts.pdf

source(here::here("reproduce", "human", "_human_helpers.R"))

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

# ---- 1. per-sample fluoroquinolone exposure --------------------------------
# Data_S4 lists every medication a patient received in the two days before each
# stool sample, one row per drug/route, pre-classified into the four categories
# this study uses. A sample counts as fluoroquinolone-exposed if any row in its
# window is in the fluoroquinolones class (ciprofloxacin or levofloxacin, oral or
# IV). This is deliberately NOT the mutually exclusive hierarchy used for the E2b/c
# histograms: here broad-spectrum and fluoroquinolone exposure are two separate
# terms in one model, so a sample can carry both.
med <- read_csv(released("Data_S4_Medication_Exposures_in_the_Two_Days_Prior_to_Stool_Sample_Collection.csv"),
                show_col_types = FALSE)

sample_exposures <- med |>
  group_by(sampleid) |>
  summarise(fluoroquinolone = any(drug_category_for_this_study == "fluoroquinolones"),
            broad_spectrum  = any(drug_category_for_this_study == "broad_spectrum"),
            .groups = "drop")

# ---- 2. model data ---------------------------------------------------------
# Same cohort and same factor coding as 10_fit_diversity_models.R. No /100 rescale
# is needed here because no dietary (gram-scale) predictor enters this model.
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) |>
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid)) |>
  inner_join(sample_exposures, by = "sampleid")

# Sanity check: the released `empirical` flag must be the broad-spectrum exposure
# re-derived from Data_S4. If this ever fails the two antibiotic terms are not
# measuring what the panel says they measure, so stop rather than fit.
stopifnot(identical(meta$empirical, meta$broad_spectrum))
message(sprintf("cohort: %d samples, %d patients | broad-spectrum %d, fluoroquinolone %d, both %d",
                nrow(meta), n_distinct(meta$pid),
                sum(meta$empirical), sum(meta$fluoroquinolone),
                sum(meta$empirical & meta$fluoroquinolone)))

# ---- 3. fit the simplified model -------------------------------------------
# The F2 prior block, copied verbatim from 10 (as 13d/16 also do): a general N(0,1)
# on every coefficient, tight N(0, 0.1) on the two nutrition-support flags, N(0, 0.5)
# on broad-spectrum exposure, and informative N(2, 0.1) intercepts for the three
# conditioning intensities. The one new coefficient here, fluoroquinolone exposure,
# is left on the general N(0,1): giving it the broad-spectrum term's tighter 0.5
# would be an invention, and the fit is unchanged either way (checked, see ledger).
diversity_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
}

# `0 + intensity` gives one intercept per conditioning intensity (the informative
# N(2, 0.1) priors in diversity_priors() are written for those cell means), and the
# two varying intercepts absorb between-patient differences and the strong
# time-since-transplant signal that panel e then reads back out.
fit <- brm(
  bf(log(simpson_reciprocal) ~ 0 + intensity + empirical + fluoroquinolone +
       TPN + EN + (1 | pid) + (1 | timebin)),
  data = meta, prior = diversity_priors(),
  warmup = 1000, iter = 3000, chains = 4, cores = 4,
  seed = 123, silent = 2, refresh = 0,
  control = list(adapt_delta = 0.99), backend = brms_backend,
  file = cache_path("E2de_fit_fluoroquinolone"), file_refit = "on_change")

# Coefficient table via fixef() (posterior mean + 2.5/97.5 percentiles), matching
# every other forest in the paper and the E2 legend's error-bar definition.
results_df <- fixef(fit, probs = c(0.025, 0.975)) |>
  as.data.frame() |>
  rownames_to_column("term") |>
  transmute(effect = "fixed", term, estimate = Estimate,
            conf.low = Q2.5, conf.high = Q97.5,
            n_samples = nrow(meta), n_patients = n_distinct(meta$pid))
write_csv(results_df, cache_path("E2de_results_df_fluoroquinolone.csv"))

# ---- 4. E2d: forest of the four exposure coefficients ----------------------
# The three intensity cell means (~2 on the ln-diversity scale) are intercepts, not
# effects, so they are dropped from the panel as in the published figure.
term_labels <- c(empiricalTRUE       = "broad-spectrum antibiotics",
                 fluoroquinoloneTRUE = "prophylactic fluoroquinolones",
                 TPNTRUE             = "total parenteral nutrition (TPN)",
                 ENTRUE              = "enteral nutrition (EN)")

forest_df <- results_df |>
  filter(term %in% names(term_labels)) |>
  mutate(label = factor(term_labels[term], levels = rev(term_labels)),
         # red where the 95% credible interval clears zero
         significant = (conf.low > 0) | (conf.high < 0))

# The published panel prints the model formula above the forest; keep it, wrapped
# the same way, hanging-indented under the outcome.
formula_caption <- paste(
  "ln(alpha-diversity) ~ conditioning +",
  "                 total parenteral nutrition +",
  "                 enteral nutrition +",
  "                 broad-spectrum antibiotics +",
  "                 prophylactic fluoroquinolones +",
  "                 (1|patient_id) +",
  "                 (1|week)", sep = "\n")

e2d <- ggplot(forest_df, aes(estimate, label, color = significant)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  linewidth = 1.6, size = 0.55, fatten = 2.6) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
  # ticks every 0.2 out to +0.2, as printed, so the axis does not stop at the
  # widest interval
  scale_x_continuous(breaks = seq(-0.4, 0.2, 0.2)) +
  expand_limits(x = c(-0.4, 0.2)) +
  labs(x = "ln(diversity) change", y = NULL, title = formula_caption) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(color = "black"),
        plot.title = element_text(size = 9.5, hjust = 0, face = "plain",
                                  margin = margin(b = 10)))

# ---- 5. E2e: the week varying intercepts -----------------------------------
# ranef() returns the per-level deviations with the same mean + 95% interval
# summary as fixef(). timebin levels are 7-day bins of transplant day, so the
# bin's lower edge / 7 is the week; the earliest bin is open-ended below (it holds
# the handful of pre-hospitalisation samples), hence the "(-2]" tick label.
week_df <- ranef(fit, probs = c(0.025, 0.975))$timebin[, , "Intercept"] |>
  as.data.frame() |>
  rownames_to_column("timebin") |>
  transmute(timebin, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5,
            week = as.numeric(str_extract(timebin, "(?<=\\[)-?\\d+")) / 7) |>
  arrange(week) |>
  mutate(week_label = if_else(week == min(week), paste0("(", week, "]"), as.character(week)),
         week_label = factor(week_label, levels = week_label),
         significant = (conf.low > 0) | (conf.high < 0))

e2e <- ggplot(week_df, aes(week_label, estimate, color = significant)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  linewidth = 1.6, size = 0.55, fatten = 2.6) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
  labs(x = "Week relative to transplant", y = "Varying intercept coefficient") +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "black"))

save_panel(e2d, "E2d_fluoroquinolone_forest.pdf", width = 130, height = 90)
save_panel(e2e, "E2e_week_intercepts.pdf", width = 120, height = 80)

# ---- 6. console summary for verification against the published panels -------
message("\nE2d coefficients:")
print(forest_df |> transmute(label, estimate = round(estimate, 3),
                             conf.low = round(conf.low, 3),
                             conf.high = round(conf.high, 3), significant))
message("E2e week intercepts:")
print(week_df |> transmute(week_label, estimate = round(estimate, 3),
                           conf.low = round(conf.low, 3),
                           conf.high = round(conf.high, 3), significant))
