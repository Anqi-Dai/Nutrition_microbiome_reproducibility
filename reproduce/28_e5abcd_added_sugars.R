# Extended Fig. E5 a-d: splitting dietary "Sugars" into ADDED vs OTHER sugars and
# refitting the macronutrient diversity model.
#
# Added sugars per food come from the USDA Food Patterns Equivalents Database (FPED),
# which reports added sugars in TEASPOON EQUIVALENTS per 100 g of food. One teaspoon
# equivalent is taken as exactly 4.2 g of (total) sugars, so
#   added_sugars_g = ADD_SUGARS(tsp/100g) * 4.2 * total_weight/100   (per food item)
#   other_sugars_g = total Sugars_g - added_sugars_g                 (clamped at 0)
# The lookup prefers FPED 2015-16; the two salad-dressing codes that exist only in
# the 2017-20 release (83112400, 83113500) are filled from FPED_1720. Each per-sample
# exposure is the prior-2-day average (sum over [sdrt-2, sdrt-1] / 2), per 100 g.
#
#   E5a  added vs total sugars per sample (Pearson; n = 1009)
#   E5b  diversity model with added + other sugars (+ Fat + Fiber), each x abx
#   E5c  diversity model with added sugars only (+ Fat + Fiber)
#   E5d  diversity model with other sugars only (+ Fat + Fiber)
# Forest styling matches the F2d/E4 diversity forests.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages(library(readxl))
if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

TSP_TO_G <- 4.2

# ---- 1. added-sugars lookup (tsp eq per 100g) ------------------------------
read_add <- function(f) {
  d <- read_excel(released(f), sheet = 1)
  col <- grep("ADD_SUGARS", names(d), value = TRUE)[1]
  tibble(Food_code = as.numeric(d$FOODCODE), add_tsp = as.numeric(d[[col]]))
}
fped <- bind_rows(read_add("FPED_1516.xls"),
                  anti_join(read_add("FPED_1720.xls"), read_add("FPED_1516.xls"),
                            by = "Food_code")) |>
  distinct(Food_code, .keep_all = TRUE)

# ---- 2. per-food added / other sugars, then daily and prior-2-day exposures -
dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE) |>
  left_join(fped, by = "Food_code") |>
  mutate(add_tsp = replace_na(add_tsp, 0),
         added_sugars_g = add_tsp * TSP_TO_G * total_weight / 100,
         other_sugars_g = pmax(Sugars_g - added_sugars_g, 0))

daily <- dtb |> group_by(pid, fdrt) |>
  summarise(added_sugars = sum(added_sugars_g), other_sugars = sum(other_sugars_g),
            Fat = sum(Fat_g), Fiber = sum(Fibers_g),
            total_sugars = sum(Sugars_g), .groups = "drop")

meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)
expo <- meta |> select(pid, sdrt, sampleid) |>
  mutate(ws = sdrt - 2, we = sdrt - 1) |>
  left_join(daily, by = join_by(pid, ws <= fdrt, we >= fdrt)) |>
  group_by(sampleid) |>
  summarise(across(c(added_sugars, other_sugars, Fat, Fiber, total_sugars),
                   ~ sum(.x, na.rm = TRUE) / 2), .groups = "drop")

model_data <- meta |>
  select(pid, sampleid, simpson_reciprocal, empirical, intensity, EN, TPN, timebin) |>
  inner_join(expo, by = "sampleid") |>
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid),
         across(c(added_sugars, other_sugars, Fat, Fiber), ~ .x / 100))

# ---- 3. E5a scatter ---------------------------------------------------------
ct <- cor.test(model_data$total_sugars / 100, model_data$added_sugars)
e5a <- ggplot(model_data, aes(total_sugars / 100, added_sugars)) +
  geom_point(colour = "#3B4FA0", alpha = 0.5, size = 0.9) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.3, size = 4,
           label = sprintf("added vs total sugars\nPearson r = %.2f\np < 0.001\nn = %s samples",
                           ct$estimate, format(nrow(model_data), big.mark = ","))) +
  labs(x = "Total Sugars (per 100g)", y = "Added Sugars (per 100g)") +
  theme_classic() + theme(aspect.ratio = 1.3)
save_panel(e5a, "E5a_added_vs_total_sugars.pdf", width = 90, height = 90)
message(sprintf("E5a Pearson r = %.3f, n = %d", ct$estimate, nrow(model_data)))

# ---- 4. fit the three models ------------------------------------------------
diversity_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
}
fit_div <- function(vars, cache) {
  fs <- paste("log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
              paste(paste(vars, "empirical", sep = "*"), collapse = " + "),
              "+ (1 | pid) + (1 | timebin)")
  brm(bf(as.formula(fs)), data = model_data, prior = diversity_priors(),
      warmup = 1000, iter = 3000, chains = 4, cores = 4, seed = 123,
      silent = 2, refresh = 0, control = list(adapt_delta = 0.99),
      backend = brms_backend, file = cache_path(cache), file_refit = "on_change")
}
fits <- list(
  b = fit_div(c("added_sugars", "other_sugars", "Fat", "Fiber"), "E5_fit_added_other"),
  c = fit_div(c("added_sugars", "Fat", "Fiber"),                 "E5_fit_added"),
  d = fit_div(c("other_sugars", "Fat", "Fiber"),                 "E5_fit_other"))

# ---- 5. forests (F2d styling) ----------------------------------------------
forest_effects <- function(fit, level_order) {
  fixef(fit, probs = c(0.025, 0.975)) |> as.data.frame() |> rownames_to_column("term") |>
    transmute(term, estimate = Estimate, conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2)) |>
    mutate(clean_term = term |>
        str_replace("empiricalTRUE$", "abx") |>
        str_remove_all("TRUE$") |>
        str_replace("empiricalTRUE:", "abx * ") |>
        str_replace_all("_", " ")) |>
    filter(!str_detect(clean_term, "intensity")) |>
    mutate(is_significant = (conf.low * conf.high) > 0,
           clean_term = factor(clean_term, levels = level_order))
}
forest_plot <- function(cleaned, title) {
  shading <- cleaned |> mutate(y = as.numeric(clean_term)) |> filter(str_detect(clean_term, "\\*"))
  ggplot(cleaned, aes(estimate, clean_term)) +
    geom_rect(data = shading, aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
              fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, color = "blue", linewidth = 0.8) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                    size = 0.25, linewidth = 1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
    scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                     str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
    labs(x = "ln(diversity) change", y = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_markdown(), plot.title = element_text(hjust = 0.5, size = 10),
          aspect.ratio = 1.5)
}
ord_b <- rev(c("abx", "EN", "TPN", "abx * added sugars", "added sugars",
               "abx * other sugars", "other sugars", "abx * Fat", "Fat", "abx * Fiber", "Fiber"))
ord_c <- rev(c("abx", "EN", "TPN", "abx * added sugars", "added sugars",
               "abx * Fat", "Fat", "abx * Fiber", "Fiber"))
ord_d <- rev(c("abx", "EN", "TPN", "abx * other sugars", "other sugars",
               "abx * Fat", "Fat", "abx * Fiber", "Fiber"))

save_panel(forest_plot(forest_effects(fits$b, ord_b), "model with\nadded + other sugars"),
           "E5b_added_other_sugars_forest.pdf", width = 80, height = 110)
save_panel(forest_plot(forest_effects(fits$c, ord_c), "model with\nadded sugars"),
           "E5c_added_sugars_forest.pdf", width = 80, height = 100)
save_panel(forest_plot(forest_effects(fits$d, ord_d), "model with\nother sugars"),
           "E5d_other_sugars_forest.pdf", width = 80, height = 100)
message("wrote E5a-d panels to results/")
