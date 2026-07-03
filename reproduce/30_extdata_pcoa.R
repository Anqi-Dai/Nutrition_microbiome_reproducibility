# Extended Fig. E3: genus-level Bray-Curtis PCoA of the human stool samples and the
# food-group/abx predictors of compositional state.
#
# Refactor of R30. Three panels:
#   E3a  PCoA scatter coloured by the dominating genus (>30% relative abundance),
#        Enterococcus split into E. faecium (asv_1) vs other; deterministic
#   E3b  the same PCoA coloured by alpha diversity (inverse Simpson); deterministic
#   E3c  forest of two brms models (PCoA Axis 1 and Axis 2 as outcomes) over the
#        food groups, their abx interactions and the clinical covariates
#
# The PCoA is deterministic arithmetic off the released genus relab table; the two
# axis models are the analysis on the path and cache to intermediate_data/.
#
# The "contrasting" genus colours are hard-coded below (the original pulled them
# from an unshipped genus colour key); they reproduce the published E3a legend.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({
  library(ape)
  library(vegan)
})

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

scatter_size <- 2

meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)

asv_relab_97 <- read_csv(released("45_quality_asv_relab_pident97_genus.csv"),
                         show_col_types = FALSE) |>
  inner_join(meta |> select(sampleid, simpson_reciprocal), by = "sampleid")

# Aggregate ASVs to genus, splitting Enterococcus into E. faecium (asv_1) vs other.
genus_relab <- asv_relab_97 |>
  filter(!is.na(genus)) |>
  mutate(genus_group = case_when(
    asv_key == "asv_1" & genus == "Enterococcus" ~ "Enterococcus_faecium",
    asv_key != "asv_1" & genus == "Enterococcus" ~ "Enterococcus_other",
    TRUE ~ as.character(genus))) |>
  group_by(sampleid, genus_group) |>
  summarise(genus_relative = sum(count_relative, na.rm = TRUE), .groups = "drop") |>
  rename(genus = genus_group)

# Per-sample domination status: the single most abundant genus if it clears 30%
# relative abundance, otherwise "Diverse/None".
dominant_taxa <- genus_relab |>
  group_by(sampleid) |>
  slice_max(order_by = genus_relative, n = 1, with_ties = FALSE) |>
  ungroup()

domination_metadata_partial <- dominant_taxa |>
  mutate(domination_status = if_else(genus_relative > 0.30,
                                     as.character(genus), "Diverse/None")) |>
  select(sampleid, domination_status, top_genus_relab = genus_relative)

sample_metadata_with_domination <- asv_relab_97 |>
  distinct(sampleid) |>
  left_join(domination_metadata_partial, by = "sampleid") |>
  mutate(domination_status = replace_na(domination_status, "Diverse/None"))

domination_leaderboard <- sample_metadata_with_domination |>
  count(domination_status, sort = TRUE, name = "sample_count")

# Collapse everything below the top 10 dominators to "Other" for cleaner colouring.
top_10_statuses <- domination_leaderboard |> slice_head(n = 10) |> pull(domination_status)

sample_metadata_with_domination <- sample_metadata_with_domination |>
  mutate(domination_status_limited = if_else(
    domination_status %in% top_10_statuses, as.character(domination_status), "Other"))

# Genus relab to wide matrix, Bray-Curtis, PCoA.
genus_matrix <- genus_relab |>
  pivot_wider(id_cols = sampleid, names_from = genus,
              values_from = genus_relative, values_fill = 0) |>
  column_to_rownames("sampleid")

bray_dist <- vegdist(genus_matrix, method = "bray")
pcoa_results <- pcoa(bray_dist)
variance_explained <- pcoa_results$values$Relative_eig[1:2] * 100

plot_data <- as.data.frame(pcoa_results$vectors) |>
  rownames_to_column("sampleid") |>
  left_join(sample_metadata_with_domination, by = "sampleid") |>
  filter(!is.na(domination_status_limited)) |>
  full_join(meta |> select(sampleid, simpson_reciprocal), by = "sampleid")

# E3a: PCoA coloured by dominating genus -----------------------------------------
# The eight non-special dominators take the contrasting palette below; E. faecium,
# the diverse/none class and the collapsed "Other" class get fixed colours.
contrasting_colors <- c(
  Bifidobacterium       = "#A77097",
  Blautia               = "#EC9B96",
  Lactobacillus         = "#3B51A3",
  Bacteroides           = "#519C8C",
  Erysipelatoclostridium = "#7A6920",
  Streptococcus         = "#9FB846",
  Holdemanella          = "#67621C",
  Akkermansia           = "#377EB8"
)

final_color_vector <- c(
  contrasting_colors,
  "Enterococcus_faecium" = "#129246",
  "Diverse/None"         = "#333333",
  "Other"                = "#AFB4B5"
)

plot_e3a <- ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = domination_status_limited)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = scatter_size, alpha = 0.7) +
  scale_color_manual(values = final_color_vector, name = "Domination status") +
  labs(x = sprintf("Axis 1 (%.2f%%)", variance_explained[1]),
       y = sprintf("Axis 2 (%.2f%%)", variance_explained[2])) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "bottom") +
  guides(color = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5)) +
  coord_fixed(ratio = 1)
save_panel(plot_e3a, "E3a_pcoa_domination.pdf", width = 140, height = 180)

# E3b: the same PCoA coloured by alpha diversity ---------------------------------
plot_e3b <- ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = simpson_reciprocal)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = scatter_size, alpha = 0.7, shape = 16) +
  scale_color_gradientn(name = bquote(""~alpha*"-Diversity (inverse Simpson index)"),
                        values = c(0, 0.4, 1),
                        colours = c("#1E90FF", "#FFFF00", "#FF3030"),
                        trans = "log", breaks = c(2, 5, 10, 20, 40)) +
  labs(x = sprintf("Axis 1 (%.2f%%)", variance_explained[1]),
       y = sprintf("Axis 2 (%.2f%%)", variance_explained[2])) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.35, "cm")) +
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                barwidth = 10, barheight = 0.6)) +
  coord_fixed(ratio = 1)
save_panel(plot_e3b, "E3b_pcoa_diversity.pdf", width = 140, height = 160)

# E3c: Axis 1 / Axis 2 as outcomes of food groups x abx --------------------------
meta_Axis <- meta |>
  mutate(intensity = factor(intensity, levels = c("nonablative", "reduced", "ablative")),
         pid = factor(pid)) |>
  mutate(across(starts_with("fg_"), ~ .x / 100)) |>
  inner_join(plot_data |> select(sampleid, Axis.1, Axis.2,
                                 domination_status, domination_status_limited),
             by = "sampleid") |>
  mutate(timebin = cut_width(sdrt, 7, boundary = 0, closed = "left"))

all_food_vars <- meta_Axis |> select(starts_with("fg_")) |> colnames()
interaction_terms <- paste(all_food_vars, "empirical", sep = "*")

priors <-
  prior(normal(0, 1), class = "b") +
  prior(normal(0, 0.1), class = "b", coef = "TPNTRUE") +
  prior(normal(0, 0.1), class = "b", coef = "ENTRUE") +
  prior(normal(0, 0.5), class = "b", coef = "empiricalTRUE")

outcomes_to_run <- c("Axis.1", "Axis.2")

fit_axis_model <- function(current_outcome) {
  formula_string <- paste(
    current_outcome,
    "~ 0 + intensity + empirical + TPN + EN +",
    paste(all_food_vars, collapse = " + "), "+",
    paste(interaction_terms, collapse = " + "),
    "+ (1 | pid) + (1 | timebin)")
  message("Fitting E3c model for ", current_outcome, " ...")
  brm(bf(as.formula(formula_string)),
      data = meta_Axis, prior = priors,
      warmup = 1000, iter = 3000, chains = 4, cores = 4,
      seed = 123, silent = 2, refresh = 0,
      control = list(adapt_delta = 0.99),
      backend = brms_backend,
      file = cache_path(paste0("R30_fit_", make.names(current_outcome))),
      file_refit = "on_change")
}

model_fits_list <- set_names(map(outcomes_to_run, fit_axis_model), outcomes_to_run)

all_results <- map_dfr(model_fits_list, ~ tidy(.x, conf.int = TRUE) |>
                         mutate(conf.low = round(conf.low, 2),
                                conf.high = round(conf.high, 2)),
                       .id = "outcome")

# Forest, styled like the published panel.
key <- food_key()
replacement_dictionary <- setNames(key$shortname, key$fg1_name)

level_order <- rev(c(
  "abx", "EN", "TPN",
  "abx * Sweets", "Sweets", "abx * Grains", "Grains",
  "abx * Milk", "Milk", "abx * Eggs", "Eggs",
  "abx * Legumes", "Legumes", "abx * Meats", "Meats",
  "abx * Fruits", "Fruits", "abx * Oils", "Oils",
  "abx * Vegetables", "Vegetables"))

new_facet_labels <- c("Axis.1" = "PCo-1", "Axis.2" = "PCo-2")

cleaned_effects <- all_results |>
  filter(effect == "fixed") |>
  mutate(
    effect_type = if_else(str_detect(term, ":"), "Interaction", "Main Effect"),
    clean_term = term |>
      str_replace("empiricalTRUE$", "abx") |>
      str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") |>
      str_replace_all(replacement_dictionary) |>
      str_replace("empiricalTRUE:", "abx * ") |>
      str_replace_all("_", " ")) |>
  filter(!str_detect(clean_term, "intensity")) |>
  mutate(is_significant = (conf.low * conf.high) > 0,
         outcome_label = recode(outcome, !!!new_facet_labels)) |>
  mutate(clean_term = factor(clean_term, levels = level_order))

shading_df <- cleaned_effects |>
  mutate(y_numeric = as.numeric(clean_term)) |>
  filter(str_detect(clean_term, "\\*"))

plot_e3c <- ggplot(cleaned_effects, aes(x = estimate, y = clean_term)) +
  geom_rect(data = shading_df,
            aes(ymin = y_numeric - 0.5, ymax = y_numeric + 0.5, xmin = -Inf, xmax = Inf),
            fill = "#FBEADC", alpha = 0.7, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "solid", color = "blue", linewidth = 0.8) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = is_significant),
                  size = 0.25, linewidth = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_y_discrete(labels = function(x) ifelse(str_detect(x, "\\*"),
                   str_glue("<b style='color:royalblue'>{x}</b>"), as.character(x))) +
  facet_wrap(~ outcome_label) +
  labs(x = "Predicted Shift in Axis Score",
       title = "Predictors of Microbiome Compositional State", y = "") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.y = element_markdown())
save_panel(plot_e3c, "E3c_axis_effects.pdf", width = 180, height = 140)

message("E3 a,b,c panels written to results/.")
