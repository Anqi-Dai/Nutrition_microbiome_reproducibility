# Extended Data E1b: relative contribution of covariates to microbiome variation
# (RESTRICTED). Ported from R37_covariates_contribution.Rmd and its helper
# identify_covariates.R.
#
# Each covariate's "variance explained" is the vegan::envfit r^2 of that covariate
# fit against a Bray-Curtis ordination of the per-sample ASV counts. Covariates that
# are FDR-significant (BH, adjusted p < 0.05) are drawn as a horizontal bar chart,
# coloured by category. E1b is the "zoomed-in" view that excludes the dominating
# per-patient factor (pid) so the remaining effects are comparable.
#
# Two changes from R37 per request:
#   - the clinical covariates (source, intensity, age, sex, disease.simple) come from
#     the cleaned restricted df_main_clinical_outcome.rds instead of
#     R02_cleaned_clinical_outcome.rds (a drop-in: same columns);
#   - the ASV counts come from released 63_asv_count_relab_res.csv (asv_key,
#     sampleid, count, count_relative) instead of R25_asv_counts.csv (same schema).
#
# df_main is restricted, so this panel skips cleanly when restricted_data/ is absent.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({
  library(vegan)
  library(forcats)
})

df_file <- "df_main_clinical_outcome.rds"
if (!has_restricted(df_file)) {
  message("E1b skipped: restricted df_main not found (", restricted(df_file), ").")
  message("This panel needs the clinical covariates (source/intensity/age/sex/disease); ",
          "place the cleaned df_main_clinical_outcome.rds in restricted_data/.")
  quit(save = "no", status = 0)
}

results_dir <- here::here("results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

set.seed(123)  # envfit permutations + the correlated-pair tie-break are stochastic

# --- envfit helpers (inlined from identify_covariates.R) --------------------
# Drop one of any near-perfectly correlated numeric pair so envfit vectors are not
# collinear (with these covariates -- age + nine food groups -- nothing crosses the
# threshold, so this is a no-op, but it is kept to match the reference).
remove_highly_correlated <- function(df, threshold = 0.95) {
  numeric_df <- df |> select(where(is.numeric))
  if (ncol(numeric_df) < 2) return(df)
  cm <- abs(cor(numeric_df, use = "pairwise.complete.obs"))
  cm <- cm[rowSums(is.na(cm)) != ncol(cm), colSums(is.na(cm)) != nrow(cm), drop = FALSE]
  drop <- character(0)
  for (i in seq_len(ncol(cm) - 1)) for (j in (i + 1):ncol(cm)) {
    v1 <- colnames(cm)[i]; v2 <- colnames(cm)[j]
    if (!(v1 %in% drop) && !(v2 %in% drop) && cm[v1, v2] > threshold) {
      drop <- c(drop, if (sample(c(TRUE, FALSE), 1)) v1 else v2)
    }
  }
  if (length(drop)) { message("Removing correlated: ", paste(drop, collapse = ", ")) }
  select(df, -any_of(drop))
}

# envfit each covariate against the distance matrix; return r^2 (effect size) and the
# BH-adjusted permutation p-value. Rows of cov_df are samples (rownames = sampleid).
find_covariates <- function(distance_mat, cov_df, perms = 10000, threshold = 0.05) {
  cov_df <- na.omit(cov_df)
  keep <- rownames(distance_mat) %in% rownames(cov_df)
  distance_mat <- distance_mat[keep, keep]
  stopifnot(nrow(distance_mat) == nrow(cov_df))
  minimal <- remove_highly_correlated(cov_df)
  cov <- vegan::envfit(ord = as.data.frame(distance_mat), env = minimal, permutations = perms)
  tibble(
    variable_name = c(names(cov$vectors$r), names(cov$factors$r)),
    effect_size   = c(cov$vectors$r, cov$factors$r),
    pval          = c(cov$vectors$pvals, cov$factors$pvals)
  ) |>
    mutate(variable_name = str_replace_all(variable_name, "_", " "),
           adjusted_pval = p.adjust(pval, method = "BH"),
           significant   = adjusted_pval < threshold)
}

# --- Data -------------------------------------------------------------------
# Per-patient clinical covariates from the cleaned restricted df_main.
clinical_cov <- read_rds(restricted(df_file)) |>
  select(source, intensity, age, sex, disease.simple, pid)

# Per-sample covariates from released metadata, joined to the clinical set and
# collapsed to disease lineage (myeloid vs lymphoid); samples with any other disease
# are dropped, as in R37.
metadata <- read.csv(released("153_combined_META.csv")) |>
  select(sampleid, empirical, timebin, fg_egg, fg_fruit, fg_grain, fg_legume,
         fg_meat, fg_milk, fg_oils, fg_sweets, fg_veggie, TPN, EN, pid) |>
  left_join(clinical_cov, by = "pid") |>
  mutate(disease_lineage = case_when(
    disease.simple %in% c("AML", "MDS/MPN", "CML") ~ "Myeloid",
    disease.simple %in% c("NHL", "ALL", "Myeloma", "CLL", "Hodgkins") ~ "Lymphoid",
    TRUE ~ NA_character_)) |>
  filter(!is.na(disease_lineage)) |>
  select(-disease.simple)

message("E1b samples removed for non myeloid/lymphoid disease: ",
        n_distinct(read.csv(released("153_combined_META.csv"))$sampleid) - n_distinct(metadata$sampleid))

# ASV count matrix (asv rows x sample cols), Bray-Curtis distance across samples.
otu <- read.csv(released("63_asv_count_relab_res.csv")) |>
  select(-count_relative) |>
  filter(sampleid %in% metadata$sampleid) |>
  pivot_wider(values_from = count, names_from = sampleid) |>
  column_to_rownames("asv_key")
otu[is.na(otu)] <- 0
dist_mat <- as.matrix(vegdist(t(otu), method = "bray"))

# Align metadata rows to the distance matrix, sampleid to rownames.
metadata <- metadata[match(rownames(dist_mat), metadata$sampleid), ] |>
  remove_rownames() |>
  column_to_rownames("sampleid")

message("E1b cohort: ", nrow(dist_mat), " samples, ", n_distinct(metadata$pid), " patients.")

# --- envfit -----------------------------------------------------------------
covariates <- find_covariates(dist_mat, metadata, perms = 10000, threshold = 0.05)

# --- Plot (E1b: significant covariates, pid excluded) -----------------------
display_names <- c(
  timebin = "Week Relative to Transplant", empirical = "ABX exposure (last 2 days)",
  source = "Graft Source", intensity = "Conditioning Intensity",
  EN = "Enteral Nutrition", TPN = "Total Parenteral Nutrition",
  sex = "Sex", age = "Age",
  `fg egg` = "Food Group: eggs", `fg fruit` = "Food Group: fruit",
  `fg grain` = "Food Group: grain", `fg legume` = "Food Group: legume",
  `fg meat` = "Food Group: meat", `fg milk` = "Food Group: milk",
  `fg oils` = "Food Group: oils", `fg sweets` = "Food Group: sweets",
  `fg veggie` = "Food Group: veggie")

cat_levels <- c("Per-Sample Events", "Dietary Factors", "Per-Patient Factors")
color_palette <- c("Per-Sample Events" = "#C77CFF",   # purple
                   "Dietary Factors"   = "#F8766D",   # salmon
                   "Per-Patient Factors" = "#00BFC4") # teal

plot_data <- covariates |>
  filter(significant, variable_name != "pid") |>
  mutate(
    category = case_when(
      str_starts(variable_name, "fg") ~ "Dietary Factors",
      variable_name %in% c("source", "intensity", "sex", "age") ~ "Per-Patient Factors",
      variable_name %in% c("empirical", "timebin", "EN") ~ "Per-Sample Events",
      TRUE ~ "Other Clinical"),
    category = factor(category, levels = cat_levels),
    variable_name = recode(variable_name, !!!display_names))

message("E1b significant covariates (variance explained):")
print(plot_data |> arrange(desc(effect_size)) |> select(variable_name, effect_size, category))

p <- plot_data |>
  ggplot(aes(x = effect_size, y = fct_reorder(variable_name, effect_size), fill = category)) +
  geom_col() +
  scale_fill_manual(values = color_palette, breaks = cat_levels) +
  scale_x_continuous(limits = c(0, 0.08), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Relative contribution of covariates",
       x = "Microbiome Variance Explained", y = NULL, fill = "Covariate Category") +
  theme_classic(base_size = 13) +
  theme(panel.grid.major.x = element_line(colour = "grey90"),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.7, 0.35))

ggsave(file.path(results_dir, "E1b_covariates_contribution.pdf"), p,
       width = 8, height = 5, device = "pdf")
message("Wrote results/E1b_covariates_contribution.pdf")
