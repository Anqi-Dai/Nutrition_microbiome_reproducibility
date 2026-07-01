# Sensitivity analyses of the F2 food-group diversity model (Extended Data E5).
#
# Refactor of R60: six variants of the diversity model, each refit on a different
# subset or with an extra covariate, each drawn as a forest panel, then assembled.
# One shared prior block, one fitting function and one forest function carry all
# six. Reads the expanded metadata R59_meta_expanded.csv (adds disease_lineage,
# ci_cleaned_numeric, PCA, exposure_type to the F2 columns). Fits cache via brms
# `file =` under intermediate_data/.
#
# Panel map (the ggsave names from R60):
#   E5e (left)  add a malignant-disease term            -> plot_disease
#   E5e (middle)add an HCT-CI comorbidity term          -> plot_CI
#   E5e (right) add a PCA (patient-controlled analgesia) term -> plot_pca
#   E5f         remove EN/TPN-exposed patients          -> plot_forest
#   E5g         empiric-exposed samples only            -> plot_empiric
#   E5h         pre-transplant (sdrt < 0) samples only  -> plot_pre

source(here::here("reproduce", "human", "_human_helpers.R"))

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

meta <- read_csv(released("R59_meta_expanded.csv"), show_col_types = FALSE) %>%
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid)) %>%
  mutate(across(starts_with("fg_"), ~ .x / 100))

key <- food_key()
repl <- setNames(key$shortname, key$fg1_name)
fg_vars <- meta %>% select(starts_with("fg")) %>% colnames()
fg_levels <- c("Sweets", "Grains", "Milk", "Eggs", "Legumes",
               "Meats", "Fruits", "Oils", "Vegetables")

# Shared prior block; drop the EN/TPN priors when those terms are absent, and
# point the antibiotic prior at whichever coefficient names the exposure.
sens_priors <- function(abx_coef = "empiricalTRUE", drop_en_tpn = FALSE) {
  p <- prior(normal(0, 1), class = "b") +
    prior_string("normal(0, 0.5)", class = "b", coef = abx_coef) +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
  if (!drop_en_tpn) {
    p <- p + prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
      prior(normal(0, 0.1), class = "b", coef = "ENTRUE")
  }
  p
}

# Fit one variant: base diversity model with optional extra fixed term, optional
# EN/TPN, and a configurable exposure variable (empirical vs exposure_type).
fit_sensitivity <- function(data, cache_file, abx_var = "empirical",
                            abx_coef = "empiricalTRUE", extra = NULL, drop_en_tpn = FALSE) {
  en_tpn <- if (drop_en_tpn) NULL else c("TPN", "EN")
  interactions <- paste(paste(fg_vars, abx_var, sep = "*"), collapse = " + ")
  rhs <- paste(c("0", "intensity", abx_var, en_tpn, extra, interactions,
                 "(1 | pid)", "(1 | timebin)"), collapse = " + ")
  brm(bf(as.formula(paste("log(simpson_reciprocal) ~", rhs))),
      data = data, prior = sens_priors(abx_coef, drop_en_tpn),
      warmup = 1000, iter = 3000, chains = 4, cores = 4,
      seed = 123, silent = 2, refresh = 0, control = list(adapt_delta = 0.99),
      backend = brms_backend, file = cache_file, file_refit = "on_change")
}

# Fixed-effect table from fixef() (robust to the underscore food-group names).
coef_table <- function(fit) {
  fixef(fit, probs = c(0.025, 0.975)) %>%
    as.data.frame() %>% rownames_to_column("term") %>%
    transmute(term, estimate = Estimate,
              conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2))
}

make_levels <- function(prefix_terms, abx = "abx") {
  inter <- as.vector(rbind(paste(abx, "*", fg_levels), fg_levels))
  rev(c(prefix_terms, inter))
}

# Forest panel shared by every variant. abx_coef/abx_label adapt the cleaning to
# the exposure variable; extra_fix patches any covariate label (e.g. HCT-CI).
sens_forest <- function(fit, level_order, subtitle, abx_coef = "empiricalTRUE",
                        abx_label = "abx", extra_fix = identity, aspect = NULL, base_size = 11) {
  cleaned <- coef_table(fit) %>%
    mutate(clean_term = term %>%
             str_replace(paste0(abx_coef, "$"), abx_label) %>%
             str_remove_all("avg_intake_|TRUE$") %>%
             str_replace_all(repl) %>%
             str_replace(paste0(abx_coef, ":"), paste0(abx_label, " * ")) %>%
             str_replace_all("_", " ") %>%
             extra_fix()) %>%
    filter(!str_detect(clean_term, "intensity")) %>%
    mutate(is_significant = (conf.low * conf.high) >= 0,
           clean_term = factor(clean_term, levels = level_order))
  shading <- cleaned %>% mutate(y = as.numeric(clean_term)) %>% filter(str_detect(clean_term, "\\*"))
  p <- ggplot(cleaned, aes(x = estimate, y = clean_term)) +
    geom_rect(data = shading, aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
              fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, color = "blue", linewidth = 0.6) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                    size = 0.2, linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                     str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
    labs(x = "ln(diversity) change", subtitle = subtitle, y = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "none", axis.text.y = element_markdown())
  if (!is.null(aspect)) p <- p + theme(aspect.ratio = aspect)
  p
}

# E5e left: add a malignant-disease term --------------------------------------
meta_disease <- meta %>% filter(!is.na(disease_lineage)) %>% mutate(disease_lineage = factor(disease_lineage))
fit_disease <- fit_sensitivity(meta_disease, cache_path("R60_fit_disease"), extra = "disease_lineage")
plot_disease <- sens_forest(fit_disease, make_levels(c("myeloid", "abx", "EN", "TPN")),
                            "model including\na term for underlying\nmalignant disease",
                            extra_fix = function(x) str_replace(x, "disease lineageMyeloid", "myeloid"))
save_panel(plot_disease, "E5e_left_disease.pdf", width = 110, height = 125)

# E5e middle: add an HCT-CI comorbidity term ----------------------------------
meta_ci <- meta %>% filter(!is.na(ci_cleaned_numeric))
fit_ci <- fit_sensitivity(meta_ci, cache_path("R60_fit_hctci"), extra = "ci_cleaned_numeric")
plot_CI <- sens_forest(fit_ci, make_levels(c("HCT-CI", "abx", "EN", "TPN")),
                       "model including\na term for\ncomorbidities",
                       extra_fix = function(x) str_replace(x, "ci cleaned numeric", "HCT-CI"))
save_panel(plot_CI, "E5e_middle_hctci.pdf", width = 110, height = 125)

# E5e right: add a PCA (patient-controlled analgesia) term ---------------------
meta_pca <- meta %>% filter(!is.na(PCA))
fit_pca <- fit_sensitivity(meta_pca, cache_path("R60_fit_pca"), extra = "PCA")
plot_pca <- sens_forest(fit_pca, make_levels(c("PCA", "abx", "EN", "TPN")),
                        "model including\na term for patient-\ncontrolled analgesia")
save_panel(plot_pca, "E5e_right_pca.pdf", width = 100, height = 100)

# E5f: exclude any patient ever exposed to EN or TPN ---------------------------
pids_en_tpn <- meta %>% group_by(pid) %>% filter(any(EN | TPN)) %>% ungroup() %>% pull(pid) %>% unique()
meta_no_en_tpn <- meta %>% filter(!pid %in% pids_en_tpn)
fit_e5f <- fit_sensitivity(meta_no_en_tpn, cache_path("R60_fit_rm_en_tpn"), drop_en_tpn = TRUE)
plot_forest <- sens_forest(fit_e5f, make_levels("abx"),
                           "excluding patients\nexposed to EN/TPN", aspect = 2.5, base_size = 10)
save_panel(plot_forest, "E5f_remove_en_tpn.pdf", width = 100, height = 100)

# E5g: empiric-exposed samples only -------------------------------------------
meta_empiric <- meta %>%
  filter(exposure_type %in% c("empiric", "no_broad_spectrum_exposure")) %>%
  mutate(exposure_type = factor(exposure_type, levels = c("no_broad_spectrum_exposure", "empiric")))
fit_empiric <- fit_sensitivity(meta_empiric, cache_path("R60_fit_empiric"),
                               abx_var = "exposure_type", abx_coef = "exposure_typeempiric")
plot_empiric <- sens_forest(fit_empiric, make_levels(c("empiric only", "EN", "TPN"), abx = "empiric only"),
                            "excluding samples\ncollected during\ntreatment of infection",
                            abx_coef = "exposure_typeempiric", abx_label = "empiric_only")
save_panel(plot_empiric, "E5g_empiric_only.pdf", width = 100, height = 100)

# E5h: pre-transplant (sdrt < 0) samples only ---------------------------------
meta_pre <- meta %>% filter(sdrt < 0)
fit_pre <- fit_sensitivity(meta_pre, cache_path("R60_fit_pretransplant"))
plot_pre <- sens_forest(fit_pre, make_levels(c("abx", "EN", "TPN")),
                        "restricting to only pre-\ntransplant samples")
save_panel(plot_pre, "E5h_pretransplant.pdf", width = 100, height = 100)

message("E5 sensitivity panels written to results/.")
