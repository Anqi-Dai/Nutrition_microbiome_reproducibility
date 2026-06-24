# Extended Fig. E5 j,k: splitting oral nutritional supplements (ONS) out of sweets.
#
# Refactor of R03. The enteral-nutrition and nutritional drink/shake/powermix food
# codes (food group "9") are reclassified into a separate group "N" (= ONS), then:
#   E5j  diet timecourse of ONS vs Sweets grams over transplant day (deterministic)
#   E5k  the food-group diversity model refit with ONS as its own group, drawn as
#        the abx-interaction forest (the abx x sweets effect survives the split)
#
# E5j is deterministic arithmetic off the diet tracker. E5k fits one brms model
# (same prior structure as the F2 food-group model, extended with the ONS group +
# its abx interaction); the fit caches to intermediate_data/ and is reused on rerun.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages(library(ggrastr))

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

# Diet tracker. Reclassify ONS food codes (enteral nutrition + nutritional
# drinks/shakes, but not the PowerBar) from group 9 to a new group "N", then total
# the dehydrated weight per patient-day-group.
dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE) %>%
  mutate(Food_code = as.character(Food_code))

g9 <- dtb %>% filter(str_detect(Food_code, "^9")) %>% distinct(Food_code, description, Meal)
ons_codes <- g9 %>%
  filter(str_detect(description, "Nutrition") | str_detect(Meal, "Enteral nutrition")) %>%
  filter(description != "Nutrition bar (PowerBar)") %>%
  pull(Food_code)

total_per_group <- dtb %>%
  mutate(Food_code = if_else(Food_code %in% ons_codes,
                             str_replace(Food_code, "^9", "N"), Food_code),
         fgrp1 = str_sub(Food_code, 1, 1)) %>%
  group_by(pid, fdrt, fgrp1) %>%
  summarise(grp_tol = sum(dehydrated_weight), .groups = "drop") %>%
  mutate(fg1_name = recode(fgrp1,
                           "1" = "fg_milk", "2" = "fg_meat", "3" = "fg_egg",
                           "4" = "fg_legume", "5" = "fg_grain", "6" = "fg_fruit",
                           "7" = "fg_veggie", "8" = "fg_oils", "9" = "fg_sweets",
                           "N" = "fg_ONS"))

# E5j: ONS vs Sweets diet timecourse --------------------------------------------
# spread/gather so a patient-day with one group but not the other contributes a
# zero rather than a gap, then loess-smooth grams against transplant day.
e5j_data <- total_per_group %>%
  filter(fg1_name %in% c("fg_ONS", "fg_sweets")) %>%
  select(pid, fdrt, fg1_name, grp_tol) %>%
  pivot_wider(names_from = fg1_name, values_from = grp_tol, values_fill = 0) %>%
  pivot_longer(c(fg_ONS, fg_sweets), names_to = "fg1_name", values_to = "grp_tol") %>%
  mutate(fg1_name = factor(str_remove(fg1_name, "fg_"), levels = c("ONS", "sweets")))

plot_e5j <- ggplot(e5j_data, aes(x = fdrt, y = grp_tol)) +
  rasterise(geom_point(alpha = 0.2, size = 0.2, shape = 16), dpi = 300) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "#E41A1C",
              linewidth = 1, fill = "hotpink") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray40") +
  facet_wrap(~ fg1_name, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = seq(0, 50, 20)) +
  scale_y_sqrt() +
  labs(x = "Transplant day", y = "Grams", title = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        strip.background = element_rect(color = "white", fill = "gray91", linewidth = 1.5),
        strip.text.x = element_text(size = 10),
        axis.title = element_text(size = axis_title_size),
        aspect.ratio = 1)
save_panel(plot_e5j, "E5j_ons_sweets_timecourse.pdf", width = 110, height = 70)

# E5k: food-group diversity model with ONS split out ----------------------------
# Per-stool-sample diet = mean of the two days before collection (sum / 2), one
# row per food group, widened with zero fill. Interaction columns fg_*_e = fg_* x
# abx are built by hand (faithful to R03), so the coefficients read fg_X / fg_X_e.
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)

stb_pair <- meta %>% transmute(pid, sampleid, p1d = sdrt - 1, p2d = sdrt - 2)

mean_p2d_diet <- function(pid_, p1d_, p2d_) {
  total_per_group %>%
    filter(pid == pid_, fdrt %in% c(p1d_, p2d_)) %>%
    group_by(fg1_name) %>%
    summarise(ave_fg = sum(grp_tol) / 2, .groups = "drop")
}

mean_p2d_df <- stb_pair %>%
  pmap(function(pid, sampleid, p1d, p2d) mean_p2d_diet(pid, p1d, p2d)) %>%
  set_names(meta$sampleid) %>%
  bind_rows(.id = "sampleid") %>%
  pivot_wider(names_from = fg1_name, values_from = ave_fg, values_fill = 0) %>%
  inner_join(meta %>% select(-starts_with("fg")), by = "sampleid") %>%
  mutate(timebin = cut_width(sdrt, 7, boundary = 0, closed = "left"))

fg_cols <- c("fg_milk", "fg_meat", "fg_egg", "fg_legume", "fg_grain",
             "fg_fruit", "fg_veggie", "fg_oils", "fg_sweets", "fg_ONS")

meta_ons <- mean_p2d_df %>%
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid),
         abx = if_else(as.character(empirical) == "TRUE", 1, 0),
         TPN = if_else(as.character(TPN) == "TRUE", 1, 0),
         EN  = if_else(as.character(EN)  == "TRUE", 1, 0)) %>%
  mutate(across(all_of(fg_cols), ~ .x / 100)) %>%
  mutate(across(all_of(fg_cols), ~ .x * abx, .names = "{.col}_e"))

ons_formula <- paste(
  "log(simpson_reciprocal) ~ 0 + intensity + abx + TPN + EN +",
  paste(fg_cols, collapse = " + "), "+",
  paste(paste0(fg_cols, "_e"), collapse = " + "),
  "+ (1 | pid) + (1 | timebin)")

# N(0,1) on every food / interaction term, tight priors on the clinical covariates,
# informative ln(diversity) ~ 2 intercepts per intensity (same structure as F2).
ons_priors <- map(c(fg_cols, paste0(fg_cols, "_e")),
                  ~ prior_string("normal(0, 1)", class = "b", coef = .x)) %>%
  reduce(`+`) +
  prior(normal(0, 0.1), class = "b", coef = "TPN") +
  prior(normal(0, 0.1), class = "b", coef = "EN") +
  prior(normal(0, 0.5), class = "b", coef = "abx") +
  prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
  prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
  prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")

message("Fitting ONS-reclassified food-group diversity model (E5k) ...")
ons_fit <- brm(bf(as.formula(ons_formula)),
               data = meta_ons, prior = ons_priors,
               warmup = 1000, iter = 3000, chains = 4, cores = 4,
               seed = 123, silent = 2, refresh = 0,
               control = list(adapt_delta = 0.99), sample_prior = TRUE,
               backend = brms_backend,
               file = cache_path("R03_fit_ons"), file_refit = "on_change")

# Forest of the fixed effects (intensity intercepts dropped), labelled and ordered
# like the published panel; significant rows (95% CI excludes zero) in red, the
# abx interaction rows shaded.
fg_label <- c(fg_milk = "Milk", fg_meat = "Meats", fg_egg = "Eggs",
              fg_legume = "Legumes", fg_grain = "Grains", fg_fruit = "Fruits",
              fg_veggie = "Vegetables", fg_oils = "Oils", fg_sweets = "Sweets",
              fg_ONS = "ONS")
fg_order <- c("Vegetables", "abx * Vegetables", "Oils", "abx * Oils",
              "Fruits", "abx * Fruits", "Meats", "abx * Meats",
              "Legumes", "abx * Legumes", "Eggs", "abx * Eggs",
              "Milk", "abx * Milk", "Grains", "abx * Grains",
              "Sweets", "abx * Sweets", "ONS", "abx * ONS", "TPN", "EN", "abx")

e5k <- fixef(ons_fit, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% rownames_to_column("term") %>%
  transmute(term, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5) %>%
  filter(!str_detect(term, "^intensity")) %>%
  mutate(base = str_remove(term, "_e$"),
         clean_term = case_when(
           str_detect(term, "_e$") ~ paste("abx *", fg_label[base]),
           term %in% names(fg_label) ~ unname(fg_label[term]),
           TRUE ~ term),
         is_significant = (conf.low * conf.high) > 0,
         clean_term = factor(clean_term, levels = fg_order))

e5k_shading <- e5k %>% mutate(y = as.numeric(clean_term)) %>% filter(str_detect(clean_term, "\\*"))

plot_e5k <- ggplot(e5k, aes(x = estimate, y = clean_term)) +
  geom_rect(data = e5k_shading, aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
            fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, color = "blue", linewidth = 0.8) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                  size = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                   str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
  labs(x = "ln(diversity) change", y = "") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", axis.text.y = element_markdown(), aspect.ratio = 1.5)
save_panel(plot_e5k, "E5k_ons_reclassified_effects.pdf", width = 90, height = 140)

message("E5 j,k panels written to results/.")
