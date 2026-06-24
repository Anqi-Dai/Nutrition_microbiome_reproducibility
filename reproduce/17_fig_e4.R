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

message("E4 panels written to results/.")
