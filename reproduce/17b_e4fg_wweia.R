# Extended Fig. E4 f,g: re-deriving the food-group diversity model under the WWEIA
# (What We Eat In America) food-category nomenclature instead of the FNDDS major
# groups, to show the abx*Sweets diversity finding is robust to the grouping scheme.
#
# Each food code's WWEIA Category code comes from the FNDDS "At A Glance" tables
# (2015-16 + 2019-20), mapped to the 13 WWEIA major groups by code range:
#   1xxx Milk and Dairy | 2xxx Protein Foods | 3xxx Mixed Dishes | 4xxx Grains |
#   5xxx Snacks and Sweets | 6000-6399 Fruit | 6400-6999 Vegetables |
#   7000-7699 Beverages | 7700-7899 Water | 8000-8399 Fats and Oils |
#   8400-8799 Condiments and Sauces | 8800-8999 Sugars | 9xxx (+missing) Other
#
# The food-code -> WWEIA assignment is aligned with the collaborator's WWEIA pipeline
# on two points. (1) Enteral formulas (food codes 95107xxx: Suplena, Vital, Vivonex,
# Osmolite, Glucerna) are absent from the At A Glance tables; their WWEIA category is
# 7208 "Nutritional beverages" (mid-level Sweetened Beverages), so they are remapped to
# 7208 and belong to BEVERAGES -- and to Sweet Beverages in the E4g split. (2) The two
# 2019-20-only salad-dressing codes (83112400, 83113500) resolve via the 2019-20 table
# to WWEIA 8012 -> Fats and Oils.
#
#   E4f  WWEIA nomenclature (13 groups, each crossed with abx)
#   E4g  same, but Beverages split following the collaborator's "cat5" scheme.
#        Non-coffee/tea beverages go by WWEIA sub-code:
#          Sweet Beverages  = 100% juice + sweetened (7000-7099, 7200-7299, incl. 7208
#                             enteral)
#          ASB              = diet beverages (7100-7199)
#        Coffee & tea (WWEIA 7302/7304) is split PER FOOD CODE by FPED added sugars:
#          added sugars > 0                       -> Sweet Beverages
#          else diet/artificial FPED description  -> ASB
#          else                                   -> Unsweetened Coffee and Tea
#        (Water stays its own group.)
#
# With those alignments both panels reproduce the collaborator's fit.
# Model = the F2d diversity model with these groups; red = 95% CrI clear of zero.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages(library(readxl))
if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

# ---- FNDDS food code -> WWEIA category code --------------------------------
read_fndds <- function(f) {
  d <- read_excel(released(f), sheet = 1, skip = 1)
  names(d)[1:4] <- c("Food_code", "desc", "wweia_code", "wweia_desc")
  d |> transmute(Food_code = as.numeric(Food_code), wweia_code = as.integer(wweia_code))
}
fndds <- bind_rows(
  read_fndds("2015-2016 FNDDS At A Glance - FNDDS Nutrient Values.xlsx"),
  anti_join(read_fndds("2019-2020 FNDDS At A Glance - FNDDS Nutrient Values.xlsx"),
            read_fndds("2015-2016 FNDDS At A Glance - FNDDS Nutrient Values.xlsx"),
            by = "Food_code")) |> distinct(Food_code, .keep_all = TRUE)

# ---- FPED added sugars + description (for the E4g coffee/tea split) ---------
read_fped <- function(f) {
  d <- read_excel(released(f), sheet = 1)
  add_col <- grep("ADD_SUGARS", names(d), value = TRUE)[1]
  tibble(Food_code = as.numeric(d$FOODCODE), add_tsp = as.numeric(d[[add_col]]),
         fped_desc = as.character(d$DESCRIPTION))
}
fped <- bind_rows(read_fped("FPED_1516.xls"),
                  anti_join(read_fped("FPED_1720.xls"), read_fped("FPED_1516.xls"),
                            by = "Food_code")) |> distinct(Food_code, .keep_all = TRUE)
# diet / artificial-sweetener description pattern (the collaborator's regex)
asb_pattern <- regex(
  "low[- ]calorie sweet|sugar substitute|artificial(ly)? sweet(ened|ener)|sugar[- ]free|non[- ]caloric|\\bdiet\\b",
  ignore_case = TRUE)

wweia_major <- function(code) case_when(
  is.na(code) ~ "fg_other",
  code < 2000 ~ "fg_milk_dairy",   code < 3000 ~ "fg_protein",
  code < 4000 ~ "fg_mixed_dishes", code < 5000 ~ "fg_grains",
  code < 6000 ~ "fg_snacks_sweets",code < 6400 ~ "fg_fruit",
  code < 7000 ~ "fg_vegetables",   code < 7700 ~ "fg_beverages",
  code < 8000 ~ "fg_water",        code < 8400 ~ "fg_fats_oils",
  code < 8800 ~ "fg_condiments",   code < 9000 ~ "fg_sugars", TRUE ~ "fg_other")
dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE) |>
  left_join(fndds, by = "Food_code") |>
  left_join(fped, by = "Food_code") |>
  # enteral formulas -> WWEIA 7208 (Nutritional beverages -> Beverages / Sweet Beverages)
  mutate(wweia_code = if_else(
    is.na(wweia_code) & str_detect(coalesce(description, ""), regex("^enteral formula", ignore_case = TRUE)),
    7208L, wweia_code)) |>
  mutate(grp_f = wweia_major(wweia_code),
         # coffee/tea (7302/7304) sweetener class from FPED added sugars + description
         ct_class = case_when(
           !is.na(add_tsp) & add_tsp > 0                 ~ "fg_sweet_bev",
           str_detect(coalesce(fped_desc, ""), asb_pattern) ~ "fg_asb",
           TRUE                                          ~ "fg_unsweet_ct"),
         grp_g = case_when(
           grp_f != "fg_beverages"          ~ grp_f,
           wweia_code %in% c(7302L, 7304L)  ~ ct_class,          # coffee & tea, FPED-split
           wweia_code >= 7100 & wweia_code < 7200 ~ "fg_asb",    # diet beverages
           TRUE                             ~ "fg_sweet_bev"))   # juice + sweetened (incl 7208)
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)

label_dict <- c(
  fg_milk_dairy = "Milk and Dairy", fg_protein = "Protein Foods",
  fg_mixed_dishes = "Mixed Dishes", fg_grains = "Grains",
  fg_snacks_sweets = "Snacks and Sweets", fg_fruit = "Fruit",
  fg_vegetables = "Vegetables", fg_beverages = "Beverages", fg_water = "Water",
  fg_fats_oils = "Fats and Oils", fg_condiments = "Condiments & Sauces",
  fg_sugars = "Sugars", fg_other = "Other", fg_sweet_bev = "Sweet Beverages",
  fg_asb = "ASB", fg_unsweet_ct = "Unsweetened Coffee and Tea")

diversity_priors <- function() {
  prior(normal(0, 1), class = "b") +
    prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
    prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
    prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE") +
    prior(normal(2, 0.1), class = "b", coef = "intensityablative") +
    prior(normal(2, 0.1), class = "b", coef = "intensityreduced") +
    prior(normal(2, 0.1), class = "b", coef = "intensitynonablative")
}

build_and_fit <- function(grp_col, cache) {
  per_food <- dtb |> rename(fg = all_of(grp_col)) |>
    group_by(pid, fdrt, fg) |> summarise(w = sum(dehydrated_weight), .groups = "drop")
  expo <- meta |> select(pid, sdrt, sampleid) |> mutate(ws = sdrt - 2, we = sdrt - 1) |>
    left_join(per_food, by = join_by(pid, ws <= fdrt, we >= fdrt)) |> filter(!is.na(fg)) |>
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

forest <- function(fit, groups_top_to_bottom, title) {
  order <- rev(c("abx", "EN", "TPN",
    as.vector(rbind(paste("abx *", groups_top_to_bottom), groups_top_to_bottom))))
  res <- fixef(fit, probs = c(0.025, 0.975)) |> as.data.frame() |> rownames_to_column("term") |>
    transmute(term, estimate = Estimate, conf.low = round(Q2.5, 2), conf.high = round(Q97.5, 2)) |>
    mutate(clean_term = term |>
        str_replace("empiricalTRUE$", "abx") |> str_remove_all("TRUE$") |>
        str_replace_all(label_dict) |> str_replace("empiricalTRUE:", "abx * ") |>
        str_replace_all("_", " ")) |>
    filter(!str_detect(clean_term, "intensity")) |>
    mutate(is_significant = (conf.low * conf.high) > 0,
           clean_term = factor(clean_term, levels = order))
  shading <- res |> mutate(y = as.numeric(clean_term)) |> filter(str_detect(clean_term, "\\*"))
  ggplot(res, aes(estimate, clean_term)) +
    geom_rect(data = shading, aes(ymin = y - 0.5, ymax = y + 0.5, xmin = -Inf, xmax = Inf),
              fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, color = "blue", linewidth = 0.8) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                    size = 0.2, linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
    scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                     str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
    labs(x = "ln (diversity) change", y = NULL, title = title) +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_markdown(size = 8), plot.title = element_text(hjust = 0.5, size = 10),
          aspect.ratio = 2.4)
}

fit_f <- build_and_fit("grp_f", "E4_fit_wweia")
fit_g <- build_and_fit("grp_g", "E4_fit_wweia_bevsplit")

groups_f <- c("Beverages", "Condiments & Sauces", "Fats and Oils", "Fruit", "Grains",
              "Milk and Dairy", "Mixed Dishes", "Other", "Protein Foods",
              "Snacks and Sweets", "Sugars", "Vegetables", "Water")
groups_g <- c("Sweet Beverages", "ASB", "Unsweetened Coffee and Tea", "Condiments & Sauces",
              "Fats and Oils", "Fruit", "Grains", "Milk and Dairy", "Mixed Dishes", "Other",
              "Protein Foods", "Snacks and Sweets", "Sugars", "Vegetables", "Water")

save_panel(forest(fit_f, groups_f, "WWEIA nomenclature"),
           "E4f_wweia_forest.pdf", width = 95, height = 200)
save_panel(forest(fit_g, groups_g, "WWEIA nomenclature with\ncustom beverage split"),
           "E4g_wweia_bevsplit_forest.pdf", width = 100, height = 215)

sigf <- fixef(fit_f) |> as.data.frame() |> rownames_to_column("term") |>
  filter(Q2.5 * Q97.5 > 0, str_detect(term, "empirical"))
message("E4f credible (red) abx effects:"); print(sigf[, c("term", "Estimate")])
message("wrote results/E4f_wweia_forest.pdf, results/E4g_wweia_bevsplit_forest.pdf")
