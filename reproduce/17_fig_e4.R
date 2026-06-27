# Extended Fig. E4 a,b,c,d, from the cached E4 model fits and the diet tracker.
#
# Refactor of the plotting halves of R62_E4.Rmd (a,b,d) and
# R51_cal_intake_predicting_diversity.Rmd (c). Panels:
#   E4a  pairwise repeated-measures correlations of the five macronutrients
#        (rmcorr; deterministic given the tracker, computed here)
#   E4b  five-macronutrient effects forest, from R62_fit_macro5_diversity
#   E4c  caloric-intake effects forest, with vs without interaction, from the two
#        R51 caloric fits
#   E4d  top-10 sweets/beverage food codes by effective per-meal consumption
#        (dehydrated weight vs sugar content; deterministic, from the tracker)
#   E4i  late-day sweets vs macronutrient sugar, with enteral-nutrition days
#        flagged orange (deterministic, from R47_sweets_and_sugar_late.Rmd)

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages(library(rmcorr))

dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE) %>%
  mutate(Food_code = as.character(Food_code),
         starch_g = Carbohydrates_g - Fibers_g - Sugars_g)

# E4a: pairwise repeated-measures correlations ----------------------------------
macro_cols <- c("Protein", "Fat", "Starch", "Fibers", "Sugars")

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

# rmcorr over each macronutrient pair: accounts for the repeated, within-patient
# structure a plain Pearson correlation would ignore.
rmcorr_results <- combn(macro_cols, 2, simplify = FALSE) %>%
  map(function(pair) {
    res <- rmcorr(participant = pid, measure1 = pair[1], measure2 = pair[2],
                  dataset = daily_macros %>% mutate(pid = factor(pid)))
    tibble(measure1 = pair[1], measure2 = pair[2], r = res$r, p_value = res$p)
  }) %>%
  list_rbind()

plot_e4a <- rmcorr_results %>%
  mutate(measure1 = factor(measure1, levels = macro_cols),
         measure2 = factor(measure2, levels = macro_cols)) %>%
  ggplot(aes(x = measure1, y = measure2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(r, 2)), color = "black", size = 2.5) +
  scale_fill_gradient2(low = "#15294B", mid = "yellow", high = "red",
                       midpoint = 0.55, limit = c(0, 1), name = "r") +
  scale_y_discrete(limits = rev(macro_cols)) +
  labs(title = "Pairwise Repeated-Measures Correlations") +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(), legend.position = "right")
save_panel(plot_e4a, "E4a_macro_rmcorr.pdf", width = 90, height = 90)

# E4b: five-macronutrient effects forest ----------------------------------------
macro5_results <- read_csv(cache_path("R62_results_df_macro5.csv"), show_col_types = FALSE)

macro5_level_order <- rev(c(
  "abx", "EN", "TPN",
  "abx * Sugars", "Sugars", "abx * Fiber", "Fiber", "abx * Fat", "Fat",
  "abx * Starch", "Starch", "abx * Protein", "Protein"))

macro5_clean <- macro5_results %>%
  filter(effect == "fixed") %>%
  mutate(clean_term = case_when(
    term == "empiricalTRUE" ~ "abx",
    term == "ENTRUE" ~ "EN",
    term == "TPNTRUE" ~ "TPN",
    term == "empiricalTRUE:avg_intake_Starch" ~ "abx * Starch",
    term == "empiricalTRUE:avg_intake_Fat" ~ "abx * Fat",
    term == "empiricalTRUE:avg_intake_Fibers" ~ "abx * Fiber",
    term == "empiricalTRUE:avg_intake_Protein" ~ "abx * Protein",
    term == "empiricalTRUE:avg_intake_Sugars" ~ "abx * Sugars",
    term == "avg_intake_Starch" ~ "Starch",
    term == "avg_intake_Fat" ~ "Fat",
    term == "avg_intake_Fibers" ~ "Fiber",
    term == "avg_intake_Protein" ~ "Protein",
    term == "avg_intake_Sugars" ~ "Sugars",
    TRUE ~ term)) %>%
  filter(!str_detect(clean_term, "intensity")) %>%
  mutate(is_significant = (conf.low * conf.high) >= 0,
         clean_term = factor(clean_term, levels = macro5_level_order))

macro5_shading <- macro5_clean %>%
  mutate(y_numeric = as.numeric(clean_term)) %>%
  filter(str_detect(clean_term, "\\*"))

plot_e4b <- ggplot(macro5_clean, aes(x = estimate, y = clean_term)) +
  geom_rect(data = macro5_shading,
            aes(ymin = y_numeric - 0.5, ymax = y_numeric + 0.5, xmin = -Inf, xmax = Inf),
            fill = "#d8dcc8", alpha = 0.7, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "solid", color = "blue", linewidth = 0.8) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                  size = 0.25, linewidth = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                   str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
  labs(x = "ln(diversity) change", y = "",
       subtitle = "Predictor associations with\nall five macronutrients included") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.y = element_markdown())
save_panel(plot_e4b, "E4b_macro5_effects.pdf", width = 90, height = 140)

# E4c: caloric-intake effects forest, with vs without interaction ---------------
cal_int   <- read_csv(cache_path("R51_results_df_cal_interaction.csv"), show_col_types = FALSE)
cal_noint <- read_csv(cache_path("R51_results_df_cal_nointeraction.csv"), show_col_types = FALSE)

cal_clean <- bind_rows(
  cal_int %>% mutate(model_type = "Model with Interaction"),
  cal_noint %>% mutate(model_type = "Model without Interaction")) %>%
  filter(effect == "fixed") %>%
  mutate(clean_term = term %>%
           str_replace("empiricalTRUE$", "abx") %>%
           str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") %>%
           str_replace("empiricalTRUE:", "abx * ") %>%
           str_replace_all("_", " ")) %>%
  filter(!str_detect(clean_term, "intensity")) %>%
  mutate(is_significant = (conf.low * conf.high) > 0)

cal_shading <- cal_clean %>%
  mutate(y_numeric = as.numeric(as.factor(clean_term))) %>%
  filter(str_detect(clean_term, "\\*"))

plot_e4c <- ggplot(cal_clean, aes(x = estimate, y = clean_term)) +
  geom_rect(data = cal_shading,
            aes(ymin = y_numeric - 0.5, ymax = y_numeric + 0.5, xmin = -Inf, xmax = Inf),
            fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "solid", color = "blue", linewidth = 0.8) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                  size = 0.25, linewidth = 1) +
  facet_wrap(~ model_type) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                   str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
  labs(x = "ln(diversity) change", y = "",
       title = "Models with or without interaction term") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.y = element_markdown(),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold"))
save_panel(plot_e4c, "E4c_caloric_effects.pdf", width = 140, height = 90)

# E4d: top-10 sweets / beverage food codes --------------------------------------
# Sum dehydrated weight and sugar across the dataset for the top-10 most-consumed
# food group "9" items, then divide by the number of meals each appears in to get
# the effective per-meal portion.
top_sweets <- dtb %>%
  mutate(fgrp1 = str_sub(Food_code, 1, 1)) %>%
  filter(fgrp1 == "9") %>%
  group_by(Food_code, description) %>%
  summarise(total_per_code = sum(dehydrated_weight, na.rm = TRUE),
            total_sugar_per_code = sum(Sugars_g, na.rm = TRUE),
            .groups = "drop") %>%
  slice_max(order_by = total_per_code, n = 10)

meal_counts <- dtb %>%
  filter(Food_code %in% top_sweets$Food_code) %>%
  distinct(Food_code, pid, fdrt, Meal) %>%
  count(Food_code, name = "meal_n")

e4d_data <- top_sweets %>%
  inner_join(meal_counts, by = "Food_code") %>%
  mutate(Total = total_per_code / meal_n,
         `Sugar content` = total_sugar_per_code / meal_n) %>%
  select(description, Total, `Sugar content`) %>%
  pivot_longer(c(Total, `Sugar content`), names_to = "metric", values_to = "grams") %>%
  mutate(metric = factor(metric, levels = c("Sugar content", "Total")),
         description = case_when(
           str_detect(str_to_lower(description), "plus") ~ "Nutritional drink or shake,\nready-to-drink (high calorie)",
           str_detect(str_to_lower(description), "ready-to-drink") ~ "Nutritional drink or shake,\nready-to-drink",
           str_detect(str_to_lower(description), "enteral") ~ "Enteral (tube feed) formula",
           str_detect(str_to_lower(description), "sports drink") ~ "Sports drink",
           TRUE ~ str_wrap(description, width = 35)))

y_axis_order <- e4d_data %>% filter(metric == "Total") %>% arrange(grams) %>% pull(description)
e4d_data <- e4d_data %>% mutate(description = factor(description, levels = y_axis_order))

plot_e4d <- ggplot(e4d_data, aes(x = grams, y = description, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Sugar content" = "#ffbcdc", "Total" = "#db2589")) +
  labs(x = "Effective per-meal\naverage consumption\nin dehydrated weight (grams)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.75, 0.4), legend.title = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"))
save_panel(plot_e4d, "E4d_top_sweets.pdf", width = 110, height = 110)

# E4j: marginal sweets effect during antibiotics, by conditioning intensity ------
# From the cached conditional draws of the three-way model (16). Each row is the
# posterior of the sweets slope within antibiotic-exposed patients for one
# intensity; significant (CI excludes zero) rows are drawn red.
e4j_plot_data <- read_csv(cache_path("R36_e4j_conditional_draws.csv"), show_col_types = FALSE) %>%
  pivot_longer(everything(), names_to = "term_label", values_to = "estimate") %>%
  group_by(term_label) %>%
  mean_qi(estimate) %>%
  rename(conf.low = .lower, conf.high = .upper) %>%
  mutate(term_label = fct_relevel(term_label, "Nonmyeloablative", "Reduced Intensity", "Myeloablative"),
         is_significant = (conf.low * conf.high) > 0)

plot_e4j <- ggplot(e4j_plot_data, aes(x = estimate, y = term_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                 height = 0.2, linewidth = 0.8) +
  geom_point(aes(color = is_significant), size = 4) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray50"), guide = "none") +
  labs(x = "Estimated Change in ln(Diversity)", y = "") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), plot.title.position = "plot")
save_panel(plot_e4j, "E4j_sweets_by_intensity.pdf", width = 100, height = 55)

# E4i: late-day sweets vs macronutrient sugar, enteral days flagged --------------
# Refactor of R47_sweets_and_sugar_late.Rmd. Two scatter+loess panels over the
# full list of patient-days with any intake: per-day dehydrated weight of the
# "9" (sugars/sweets/beverages) food group, and per-day total macronutrient
# sugar. Days that include an enteral-nutrition meal are coloured orange; the
# point being that the late-day sweets uptick is driven by enteral formula while
# real sugar intake keeps falling.
e4i_orange <- "#D55E00"
e4i_scatter_alpha <- 0.4

all_patient_days_with_intake <- dtb %>%
  filter(dehydrated_weight > 0) %>%
  distinct(pid, fdrt)

daily_sweets_intake <- dtb %>%
  filter(str_starts(Food_code, "9")) %>%
  mutate(is_enteral = str_detect(Meal, regex("Enteral nutrition", ignore_case = TRUE))) %>%
  group_by(pid, fdrt) %>%
  summarise(total_sweets_dehydrated_weight = sum(dehydrated_weight, na.rm = TRUE),
            had_enteral = any(is_enteral), .groups = "drop")

daily_sugar_summary <- dtb %>%
  group_by(pid, fdrt) %>%
  summarise(total_sugars_g = sum(Sugars_g, na.rm = TRUE),
            had_enteral_meal = any(str_detect(Meal, regex("Enteral nutrition", ignore_case = TRUE))),
            .groups = "drop") %>%
  mutate(intake_source = if_else(had_enteral_meal, "Had Enteral Nutrition", "Only Oral Intake")) %>%
  select(pid, fdrt, intake_source, total_sugars_g)

daily_intake_final <- all_patient_days_with_intake %>%
  left_join(daily_sweets_intake, by = c("pid", "fdrt")) %>%
  left_join(daily_sugar_summary, by = c("pid", "fdrt")) %>%
  mutate(total_sweets_dehydrated_weight = replace_na(total_sweets_dehydrated_weight, 0),
         total_sugars_g = replace_na(total_sugars_g, 0),
         had_enteral = replace_na(had_enteral, FALSE),
         had_enteral = factor(had_enteral, levels = c(FALSE, TRUE),
                              labels = c("Only Oral Intake", "Had Enteral Nutrition")),
         intake_source = factor(intake_source, levels = c("Only Oral Intake", "Had Enteral Nutrition")))

e4i_theme <- theme_bw(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        panel.grid = element_blank(), aspect.ratio = 1)
e4i_cols <- c("Only Oral Intake" = "gray", "Had Enteral Nutrition" = e4i_orange)

plot_e4i_sweets <- ggplot(daily_intake_final,
                          aes(fdrt, total_sweets_dehydrated_weight, color = had_enteral)) +
  geom_point(alpha = e4i_scatter_alpha, size = 1.5) +
  geom_smooth(method = "loess", formula = "y ~ x", colour = "#E41A1C", linewidth = 1, fill = "hotpink") +
  scale_color_manual(values = e4i_cols) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray40") +
  scale_y_sqrt() +
  labs(title = "FNDDS Sweets Food Group", x = "Day Relative to Transplant",
       y = "Total Dehydrated Weight (grams)") +
  e4i_theme
save_panel(plot_e4i_sweets, "E4i_sweets_foodgroup.pdf", width = 80, height = 90)

plot_e4i_sugar <- ggplot(daily_intake_final,
                         aes(fdrt, total_sugars_g, color = intake_source)) +
  geom_point(alpha = e4i_scatter_alpha, size = 1.5) +
  geom_smooth(method = "loess", formula = "y ~ x", colour = "#E41A1C", linewidth = 1, fill = "hotpink") +
  scale_color_manual(values = e4i_cols) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "gray40") +
  scale_y_sqrt() +
  labs(title = "Macronutrient Sugar", x = "Day Relative to Transplant",
       y = "Total Sugars (grams)") +
  e4i_theme
save_panel(plot_e4i_sugar, "E4i_macronutrient_sugar.pdf", width = 80, height = 90)

# Late-day enteral share, the quantitative claim behind E4i.
e4i_late_share <- daily_intake_final %>%
  filter(fdrt > 20) %>%
  count(intake_source) %>%
  mutate(perc = round(n / sum(n) * 100, 1))
message("E4i late-day (>20) intake source share:")
print(e4i_late_share)

# E4e: FNDDS two-digit subgroup effects forest ----------------------------------
# From the cached R27 fit (16). Group 9 is split into its FNDDS two-digit
# subgroups (91/92/94/95; 93 alcoholic is absent from the cohort), each crossed
# with antibiotics, alongside the usual eight food groups. Same forest style as
# E4b: peach-shaded interaction rows, blue zero line, significant effects red,
# interaction labels royalblue and bold.
fndds_results <- read_csv(cache_path("R27_results_df_fndds.csv"), show_col_types = FALSE)

# Map each model term to its published label; the FNDDS subgroups read as
# "NN (short name)", groups 1-8 as their food-group names.
fndds_base_label <- c(
  empiricalTRUE = "abx", ENTRUE = "EN", TPNTRUE = "TPN",
  avg_intake_Formulated_beverages = "95 (formulated bev)",
  avg_intake_Nonalcoholic_beverages = "92 (non-alcoholic bev)",
  avg_intake_Sugars_and_sweets = "91 (sugars & sweets)",
  avg_intake_Water = "94 (water)",
  avg_intake_fg_grain = "Grains", avg_intake_fg_milk = "Milk",
  avg_intake_fg_egg = "Eggs", avg_intake_fg_legume = "Legumes",
  avg_intake_fg_meat = "Meats", avg_intake_fg_fruit = "Fruits",
  avg_intake_fg_oils = "Oils", avg_intake_fg_veggie = "Vegetables")

fndds_label <- function(term) {
  ifelse(str_detect(term, "^empiricalTRUE:"),
         paste("abx *", fndds_base_label[str_remove(term, "^empiricalTRUE:")]),
         fndds_base_label[term])
}

# Top-to-bottom panel order (reversed for the bottom-up y axis).
fndds_level_order <- rev(c(
  "abx", "EN", "TPN",
  "abx * 95 (formulated bev)", "95 (formulated bev)",
  "abx * 92 (non-alcoholic bev)", "92 (non-alcoholic bev)",
  "abx * 91 (sugars & sweets)", "91 (sugars & sweets)",
  "abx * 94 (water)", "94 (water)",
  "abx * Grains", "Grains", "abx * Milk", "Milk",
  "abx * Eggs", "Eggs", "abx * Legumes", "Legumes",
  "abx * Meats", "Meats", "abx * Fruits", "Fruits",
  "abx * Oils", "Oils", "abx * Vegetables", "Vegetables"))

fndds_clean <- fndds_results %>%
  filter(effect == "fixed", !str_detect(term, "intensity")) %>%
  mutate(clean_term = fndds_label(term),
         is_significant = (conf.low * conf.high) >= 0,
         clean_term = factor(clean_term, levels = fndds_level_order))

fndds_shading <- fndds_clean %>%
  mutate(y_numeric = as.numeric(clean_term)) %>%
  filter(str_detect(clean_term, "\\*"))

plot_e4e <- ggplot(fndds_clean, aes(x = estimate, y = clean_term)) +
  geom_rect(data = fndds_shading,
            aes(ymin = y_numeric - 0.5, ymax = y_numeric + 0.5, xmin = -Inf, xmax = Inf),
            fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "solid", color = "blue", linewidth = 0.8) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                  size = 0.25, linewidth = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                   str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
  labs(x = "ln(diversity) change", y = "",
       subtitle = "FNDDS two-digit nomenclature") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.y = element_markdown())
save_panel(plot_e4e, "E4e_fndds_subgroup_forest.pdf", width = 115, height = 150)

# E4e (tree): FNDDS group-9 nomenclature excerpt with per-subgroup exposure -------
# The "excerpt of FNDDS nomenclature tree" schematic that accompanies the E4e
# forest: the major group 9 root and its five two-digit subgroups, each annotated
# with the number of stool samples (and patients) that recorded any intake of it in
# the 2-day-prior window. Drawn directly with ggplot segments + ggtext labels (the
# original used ggtree; this avoids that Bioconductor dependency). Exposure counts
# come from the cached model frame (16); the subgroup descriptions and examples are
# FNDDS nomenclature annotations carried from the manuscript figure.
grey_x <- "#b3b3b3"
grey_eg <- "#8c8c8c"

fndds_exposure <- read_csv(cache_path("R27_fndds_exposure_counts.csv"), show_col_types = FALSE)

# y is top-to-bottom in display order (91 at the top); 93 carries no exposure row
# because alcoholic beverages never appear in the cohort, so it shows 0 / 0.
tree_leaves <- tribble(
  ~code, ~food_category,           ~desc,                                  ~example,                              ~y,
  "91",  "Sugars_and_sweets",      "Sugars and sweets",                    "e.g. hard candy, syrup, brown sugar", 5,
  "92",  "Nonalcoholic_beverages", "Nonalcoholic beverages",               "e.g. ginger ale, fruit juice drink",  4,
  "93",  "Alcoholic_beverages",    "Alcoholic beverages",                  "e.g. beer, wine",                     3,
  "94",  "Water",                  "Water noncarbonated",                  "e.g. bottled water, vitamin water",   2,
  "95",  "Formulated_beverages",   "Formulated nutrition beverages etc.",  "e.g. sports drink, nutritional shake", 1) %>%
  left_join(fndds_exposure, by = "food_category") %>%
  mutate(samples = replace_na(samples, 0), patients = replace_na(patients, 0),
         exposure_text = str_glue("{samples} samples from {patients} patients"),
         code_label = str_glue("<b>{code}</b><span style='color:{grey_x}'>XXXXXX</span>"),
         desc_label = str_glue("{desc}<br><span style='color:{grey_eg}'><i>{example}</i></span>"))

# Bracket geometry: a vertical spine with a short dash out to each leaf, and one
# connector from the root label in to the spine.
x_spine <- 1; x_leaf <- 1.45; x_root <- 0
root_y <- mean(range(tree_leaves$y))
tree_segments <- bind_rows(
  tibble(x = x_spine, xend = x_spine, y = min(tree_leaves$y), yend = max(tree_leaves$y)),
  tree_leaves %>% transmute(x = x_spine, xend = x_leaf, y = y, yend = y),
  tibble(x = x_root, xend = x_spine, y = root_y, yend = root_y))

plot_e4e_tree <- ggplot() +
  geom_segment(data = tree_segments, aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.5, color = "grey25") +
  ggtext::geom_richtext(
    data = tibble(x = x_root, y = root_y,
                  label = str_glue("<b>9</b><span style='color:{grey_x}'>XXXXXXXX</span><br>Sugar, Sweets,<br>and Beverages")),
    aes(x = x, y = y, label = label), hjust = 1, vjust = 0.5, size = 3.4,
    fill = NA, label.color = NA, lineheight = 1.1) +
  ggtext::geom_richtext(data = tree_leaves, aes(x = x_leaf + 0.05, y = y, label = code_label),
                        hjust = 0, size = 3.4, fill = NA, label.color = NA) +
  ggtext::geom_richtext(data = tree_leaves, aes(x = x_leaf + 1.9, y = y, label = desc_label),
                        hjust = 0, size = 3.2, fill = NA, label.color = NA, lineheight = 1.1) +
  ggtext::geom_richtext(data = tree_leaves, aes(x = x_leaf + 7.1, y = y, label = exposure_text),
                        hjust = 0, size = 3.2, fill = NA, label.color = NA) +
  labs(title = "Excerpt of FNDDS nomenclature tree") +
  scale_x_continuous(limits = c(-3, 14)) +
  scale_y_continuous(limits = c(0.4, 5.6)) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        plot.margin = margin(6, 6, 6, 6))
save_panel(plot_e4e_tree, "E4e_fndds_nomenclature_tree.pdf", width = 200, height = 70)

message("E4 panels written to results/.")
