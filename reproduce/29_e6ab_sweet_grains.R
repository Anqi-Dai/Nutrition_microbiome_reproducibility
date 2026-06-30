# Extended Fig. E6 a,b: splitting the FNDDS "Grains" food group into "Sweet Grains"
# and "Other Grains", to show the abx*Sweets diversity finding is robust to the
# food-grouping scheme.
#
# Grain food codes (FNDDS major group 5) are classified per food code:
#   E6a  Sweet Grains  = > 1.5 tsp added-sugar equivalents per 100 g (FPED database);
#                        Other Grains = <= 1.5 tsp  (38 zero-added-sugar grains -> Other)
#   E6b  Sweet Grains  = > 10 g macronutrient sugar per 100 g (FNDDS); Other otherwise
# The F2d food-group diversity model is then refit with the nine groups, but with
# Grains replaced by the two custom sub-groups (each crossed with antibiotics):
#   log(simpson) ~ 0 + intensity + abx + TPN + EN + (10 food groups)*abx
#                  + (1|pid) + (1|timebin)
#
# Expectation (from the response letter): Sweet Grains shows a similar negative
# abx interaction to abx*Sweets but is not credible at 95% (e.g. P<0 = 89% for the
# tsp cutoff); abx and abx*Sweets remain the red effects.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({ library(readxl); library(posterior) })
if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)
key <- food_key()

# ---- FPED added-sugars per food code (tsp eq per 100 g) ---------------------
read_add <- function(f) {
  d <- read_excel(released(f), sheet = 1)
  col <- grep("ADD_SUGARS", names(d), value = TRUE)[1]
  tibble(Food_code = as.numeric(d$FOODCODE), add_tsp = as.numeric(d[[col]]))
}
fped <- bind_rows(read_add("FPED_1516.xls"),
                  anti_join(read_add("FPED_1720.xls"), read_add("FPED_1516.xls"),
                            by = "Food_code")) |> distinct(Food_code, .keep_all = TRUE)

# ---- classify each GRAIN food code (FNDDS group 5) sweet vs other ----------
dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE)
grain_class <- dtb |>
  filter(substr(as.character(Food_code), 1, 1) == "5") |>
  group_by(Food_code) |>
  summarise(sugar_per_100g = median(Sugars_g / total_weight * 100, na.rm = TRUE), .groups = "drop") |>
  left_join(fped, by = "Food_code") |>
  mutate(add_tsp = replace_na(add_tsp, 0),
         sweet_fped  = add_tsp > 1.5,           # E6a cutoff
         sweet_fndds = sugar_per_100g > 10)     # E6b cutoff
message(sprintf("grain codes: %d  | sweet by FPED-tsp: %d  | sweet by FNDDS-10g: %d",
                nrow(grain_class), sum(grain_class$sweet_fped), sum(grain_class$sweet_fndds)))

fg_map <- c("1" = "fg_milk", "2" = "fg_meat", "3" = "fg_egg", "4" = "fg_legume",
            "6" = "fg_fruit", "7" = "fg_veggie", "8" = "fg_oils", "9" = "fg_sweets")
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)

# ---- build model data with Grains split, fit, and forest --------------------
fg_dict <- c(setNames(key$shortname, key$fg1_name),
             fg_sweet_grain = "Sweet Grains", fg_other_grain = "Other Grains")
fg_order <- rev(c("abx", "EN", "TPN",
  "abx * Sweets", "Sweets", "abx * Sweet Grains", "Sweet Grains",
  "abx * Other Grains", "Other Grains", "abx * Milk", "Milk", "abx * Eggs", "Eggs",
  "abx * Legumes", "Legumes", "abx * Meats", "Meats", "abx * Fruits", "Fruits",
  "abx * Oils", "Oils", "abx * Vegetables", "Vegetables"))

diversity_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
}

build_and_fit <- function(sweet_col, cache) {
  sweet_codes <- grain_class$Food_code[grain_class[[sweet_col]]]
  per_food <- dtb |>
    mutate(d1 = substr(as.character(Food_code), 1, 1),
           fg = case_when(
             d1 == "5" & Food_code %in% sweet_codes ~ "fg_sweet_grain",
             d1 == "5" ~ "fg_other_grain",
             TRUE ~ fg_map[d1])) |>
    filter(!is.na(fg)) |>
    group_by(pid, fdrt, fg) |>
    summarise(w = sum(dehydrated_weight), .groups = "drop")
  expo <- meta |> select(pid, sdrt, sampleid) |>
    mutate(ws = sdrt - 2, we = sdrt - 1) |>
    left_join(per_food, by = join_by(pid, ws <= fdrt, we >= fdrt)) |>
    filter(!is.na(fg)) |>
    group_by(sampleid, fg) |> summarise(v = sum(w) / 2, .groups = "drop") |>
    pivot_wider(names_from = fg, values_from = v, values_fill = 0)
  data <- meta |> select(pid, sampleid, simpson_reciprocal, empirical, intensity, EN, TPN, timebin) |>
    left_join(expo, by = "sampleid") |>
    mutate(across(starts_with("fg_"), ~ replace_na(.x, 0) / 100),
           intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
           pid = factor(pid))
  fg_vars <- grep("^fg_", names(data), value = TRUE)
  fs <- paste("log(simpson_reciprocal) ~ 0 + intensity + empirical + TPN + EN +",
              paste(paste(fg_vars, "empirical", sep = "*"), collapse = " + "),
              "+ (1 | pid) + (1 | timebin)")
  brm(bf(as.formula(fs)), data = data, prior = diversity_priors(),
      warmup = 1000, iter = 3000, chains = 4, cores = 4, seed = 123,
      silent = 2, refresh = 0, control = list(adapt_delta = 0.99),
      backend = brms_backend, file = cache_path(cache), file_refit = "on_change")
}

forest <- function(fit, title) {
  res <- fixef(fit, probs = c(0.025, 0.975)) |> as.data.frame() |> rownames_to_column("term") |>
    transmute(term, estimate = Estimate, conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2)) |>
    mutate(clean_term = term |>
        str_replace("empiricalTRUE$", "abx") |> str_remove_all("TRUE$") |>
        str_replace_all(fg_dict) |> str_replace("empiricalTRUE:", "abx * ") |>
        str_replace_all("_", " ")) |>
    filter(!str_detect(clean_term, "intensity")) |>
    mutate(is_significant = (conf.low * conf.high) > 0,
           clean_term = factor(clean_term, levels = fg_order))
  shading <- res |> mutate(y = as.numeric(clean_term)) |> filter(str_detect(clean_term, "\\*"))
  ggplot(res, aes(estimate, clean_term)) +
    geom_rect(data = shading, aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
              fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, color = "blue", linewidth = 0.8) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                    size = 0.22, linewidth = 0.9) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
    scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                     str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
    labs(x = "ln(diversity) change", y = NULL, title = title) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_markdown(), plot.title = element_text(hjust = 0.5, size = 10),
          aspect.ratio = 2)
}

fit_a <- build_and_fit("sweet_fped",  "E6_fit_sweetgrain_fped")
fit_b <- build_and_fit("sweet_fndds", "E6_fit_sweetgrain_fndds")

save_panel(forest(fit_a, "\"Sweet Grains\" > 1.5 tsp added sugar\nequivalents per FPED database"),
           "E6a_sweet_grains_fped.pdf", width = 85, height = 150)
save_panel(forest(fit_b, "\"Sweet Grains\" > 10 g sugar\nper FNDDS database"),
           "E6b_sweet_grains_fndds.pdf", width = 85, height = 150)

# posterior P(abx*Sweet Grains < 0), as quoted in the response letter
pneg <- function(fit) {
  d <- as_draws_df(fit)
  v <- grep("b_empiricalTRUE.fg_sweet_grain", variables(d), value = TRUE)[1]
  mean(d[[v]] < 0)
}
message(sprintf("P(abx*Sweet Grains < 0): FPED-tsp = %.0f%%, FNDDS-10g = %.0f%%",
                100 * pneg(fit_a), 100 * pneg(fit_b)))
message("wrote results/E6a_sweet_grains_fped.pdf, results/E6b_sweet_grains_fndds.pdf")
