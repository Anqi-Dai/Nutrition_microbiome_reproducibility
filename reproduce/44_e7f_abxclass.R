# Extended Fig. E7f: E. faecium (asv_1) CLR food-group effects stratified by
# antibiotic CLASS, not just binary exposure.
#
# Refactor of R46. Each stool sample's broad-spectrum exposure (prior 2 days) is
# classed as:
#   anaerobe-targeting  piperacillin-tazobactam, meropenem, metronidazole,
#                       linezolid, or ORAL vancomycin
#   anaerobe-sparing    cefepime
#   non-broad-spectrum  everything else (the reference)
# from the medication-exposure table (Data_S4 == the cleaned 191 table). A sample
# is anaerobe-targeting if it has any targeting drug, else anaerobe-sparing if any
# sparing drug, else non-broad-spectrum.
#
# The model is the F4b CLR model with the binary `empirical` replaced by this
# 3-level exposure_class_category:
#   asv_1_clr ~ 0 + intensity + exposure_class_category + TPN + EN
#               + (nine food groups) * exposure_class_category + (1|pid) + (1|timebin)
# E7f forests the food-group terms per class: the non-broad-spectrum panel is the
# main food effect (reference), the other two panels are the interaction terms
# (difference from reference). Red = 95% CrI clear of zero.

source(here::here("reproduce", "human", "_human_helpers.R"))

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)
key <- food_key()

anaerobic_targeting <- c("piperacillin_tazobactam", "meropenem", "metronidazole", "linezolid")
anaerobic_sparing   <- c("cefepime")

# ---- 1. classify each sample's broad-spectrum exposure ---------------------
med <- read_csv(released("Data_S4_Medication_Exposures_in_the_Two_Days_Prior_to_Stool_Sample_Collection.csv"),
                show_col_types = FALSE) |>
  mutate(broad_spectrum_subtype = case_when(
    drug_name_clean %in% anaerobic_targeting ~ "anaerobic_targeting",
    drug_name_clean == "vancomycin" & route_clean == "oral" ~ "anaerobic_targeting",
    drug_name_clean %in% anaerobic_sparing ~ "anaerobic_sparing",
    TRUE ~ "non_broad_spectrum"))

sample_summary <- med |>
  group_by(sampleid, sdrt, pid) |>
  summarise(exposure_class_category = case_when(
      any(broad_spectrum_subtype == "anaerobic_targeting") ~ "anaerobic_targeting",
      any(broad_spectrum_subtype == "anaerobic_sparing")   ~ "anaerobic_sparing",
      TRUE ~ "non_broad_spectrum"), .groups = "drop")
message("exposure classes:")
print(count(sample_summary, exposure_class_category))

# ---- 2. model data ---------------------------------------------------------
meta <- read_csv(released("R59_meta_expanded.csv"), show_col_types = FALSE) |>
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid),
         across(starts_with("fg_"), ~ .x / 100)) |>
  inner_join(sample_summary, by = c("sampleid", "sdrt", "pid")) |>
  mutate(exposure_class_category = factor(exposure_class_category,
            levels = c("non_broad_spectrum", "anaerobic_sparing", "anaerobic_targeting")))
message("model n: ", nrow(meta), " samples, ", n_distinct(meta$pid), " patients")

fg_vars <- meta |> select(starts_with("fg")) |> colnames()
formula_string <- paste(
  "asv_1_clr ~ 0 + intensity + exposure_class_category + TPN + EN +",
  paste(paste(fg_vars, "exposure_class_category", sep = "*"), collapse = " + "),
  "+ (1 | pid) + (1 | timebin)")

priors <- prior(normal(0, 1), class = "b") +
  prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
  prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
  prior(normal(0, 0.5), class = "b", coef = "exposure_class_categoryanaerobic_sparing") +
  prior(normal(0, 0.5), class = "b", coef = "exposure_class_categoryanaerobic_targeting")

fit <- brm(bf(as.formula(formula_string)), data = meta, prior = priors,
           warmup = 1000, iter = 3000, chains = 4, cores = 4, seed = 123,
           silent = 2, refresh = 0, control = list(adapt_delta = 0.99),
           backend = brms_backend,
           file = cache_path("R46_fit_asv1_abxclass"), file_refit = "on_change")

results_df <- fixef(fit, probs = c(0.025, 0.975)) |> as.data.frame() |>
  rownames_to_column("term") |>
  transmute(effect = "fixed", term, estimate = Estimate,
            conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2))
write_csv(results_df, cache_path("R46_results_df_abxclass.csv"))

# ---- 3. E7f forest ---------------------------------------------------------
replacement_dictionary <- setNames(key$shortname, key$fg1_name)
plot_df <- results_df |>
  filter(effect == "fixed", str_detect(term, "fg_")) |>
  mutate(
    exposure = case_when(
      str_detect(term, "anaerobic_targeting") ~ "anaerobe-targeting\nbroad-spectrum",
      str_detect(term, "anaerobic_sparing")   ~ "anaerobe-sparing\nbroad-spectrum",
      TRUE                                     ~ "broad-spectrum\nnon-exposed"),
    food = str_extract(term, "fg_[a-z]+") |> str_replace_all(replacement_dictionary),
    significant = (conf.low > 0) | (conf.high < 0),
    exposure = factor(exposure, levels = c(
      "broad-spectrum\nnon-exposed", "anaerobe-sparing\nbroad-spectrum",
      "anaerobe-targeting\nbroad-spectrum")),
    food = factor(food, levels = rev(c(
      "Sweets", "Grains", "Milk", "Eggs", "Legumes",
      "Meats", "Fruits", "Oils", "Vegetables"))))

e7f <- ggplot(plot_df, aes(estimate, food, color = significant)) +
  geom_vline(xintercept = 0, color = "steelblue", linewidth = 0.7) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.3) +
  facet_wrap(~ exposure, nrow = 1) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
  labs(x = "CLR *E. faecium* Change per 100g Food Intake", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        axis.title.x = ggtext::element_markdown())

save_panel(e7f, "E7f_efaecium_abxclass_forest.pdf", width = 200, height = 90)
sig <- plot_df |> filter(significant) |>
  transmute(panel = str_replace_all(as.character(exposure), "\n", " "), food, estimate, conf.low, conf.high)
message("\ncredible (red) food-group effects:")
print(sig)
