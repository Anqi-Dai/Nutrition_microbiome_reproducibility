# Shared helpers for the human diversity-model figures (F2 and its diagnostics).
#
# The Bayesian model family is the contribution, so the full model code is on the
# path: 10 fits and caches the brms models, 11 draws F2 c-h from the cached fits,
# 12 draws the prior/posterior predictive diagnostics. Data root, caching, the
# food-group key and the antibiotics palette live here so the three scripts stay
# thin and consistent.

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(tidybayes)
  library(broom.mixed)
  library(posterior)
  library(ggtext)
  library(ggpubr)
  library(wesanderson)
  library(ggsci)
  library(cowplot)
  library(here)
})

# brms backend. rstan 2.32 fails to compile against the current Apple clang /
# StanHeaders, so prefer the cmdstanr backend (CmdStan compiles with its own
# toolchain config) when available, falling back to rstan otherwise.
brms_backend <- "rstan"
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  ver <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL)
  if (is.null(ver)) {
    try(cmdstanr::set_cmdstan_path(path.expand("~/.cmdstan/cmdstan-2.38.0")), silent = TRUE)
    ver <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL)
  }
  if (!is.null(ver)) brms_backend <- "cmdstanr"
}

# One data root, resolved from one environment variable (same contract as mouse).
nutrition_data_root <- function() {
  Sys.getenv("NUTRITION_DATA", unset = here::here("released_data"))
}
released <- function(file) file.path(nutrition_data_root(), file)

# Restricted tier: PHI-free but not cleared for public release, so it is
# gitignored and ships to nobody. Scripts that need it resolve from RESTRICTED_DATA
# (falling back to restricted_data/) and guard on has_restricted() so they skip
# cleanly for anyone who only has released_data/.
restricted_data_root <- function() {
  Sys.getenv("RESTRICTED_DATA", unset = here::here("restricted_data"))
}
restricted <- function(file) file.path(restricted_data_root(), file)
has_restricted <- function(file = NULL) {
  if (is.null(file)) dir.exists(restricted_data_root()) else file.exists(restricted(file))
}

# Expensive brms fits and the derived coefficient tables cache here; reused on rerun.
intermediate_dir <- function() here::here("intermediate_data")
cache_path <- function(...) file.path(intermediate_dir(), ...)

# Panel-level outputs with manuscript-matching names land in results/. The human
# panels are sized in millimetres to match the original ggsave calls in 172.
save_panel <- function(plot, file, width = 90, height = 90, units = "mm") {
  out_dir <- here::here("results")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  path <- file.path(out_dir, file)
  ggsave(path, plot = plot, width = width, height = height, units = units, device = "pdf")
  message("wrote ", path)
  invisible(path)
}

# Food-group label/colour key (fg1_name -> shortname, plus per-group hex colours).
food_key <- function() {
  read_csv(released("food_group_color_key_final.csv"), col_types = "ccccc")
}

# Antibiotics palette: Wes Anderson Royal1 first two (not exposed, exposed).
abx_palette <- wes_palette("Royal1", 2)

axis_text_size <- 10
axis_title_size <- 10

# The genus set for the whole taxon family. A genus is kept when it exceeds 0.2%
# relative abundance in more than 10% of the cohort's stool samples, which selects
# 33. The E7a heatmap models (40), the E7b correlation bar and the F4a summary all
# start from this one set, so a genus cannot appear in one panel and be missing
# from another.
prevalent_genera <- function(genus_long, n_samples) {
  genus_long |>
    filter(!is.na(genus), relab > 0.002) |>
    count(genus, name = "n") |>
    mutate(perc = round(n / n_samples * 100, 0)) |>
    filter(perc > 10) |>
    pull(genus) |>
    sort()
}

# Genus-abundance vs alpha-diversity Spearman correlations (F4a / E7b source,
# ported from 178_new_F4__code_for_Figure_4.Rmd). The original read a pre-built
# genus-count table (022_ALL173_stool_samples_genus_counts.csv); that table does
# not ship, so the genus relative abundance is rebuilt here from the released
# per-ASV genus relab (45_quality_asv_relab_pident97_genus.csv): drop the
# unassigned (NA) genus, sum count_relative to genus level, then zero-fill the
# sample x genus grid (spread/gather) exactly as the original did. Each genus
# relab is correlated against inverse-Simpson diversity, then Benjamini-Hochberg
# adjusted across the prevalent genera. A tiny seeded jitter breaks relab ties so
# cor.test can attempt exact p-values, matching 178; it is drawn over the full
# sample x genus grid before the prevalence filter, because the size of the draw
# shifts the tie-breaking and so would shift every rho.
genus_diversity_spearman <- function() {
  set.seed(1)
  meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)

  g_relab <- read_csv(released("45_quality_asv_relab_pident97_genus.csv"),
                      show_col_types = FALSE) %>%
    filter(!is.na(genus)) %>%
    group_by(sampleid, genus) %>%
    summarize(relab = sum(count_relative, na.rm = TRUE), .groups = "drop") %>%
    filter(sampleid %in% meta$sampleid) %>%
    spread("genus", "relab", fill = 0) %>%
    gather("genus", "relab", -sampleid) %>%
    inner_join(meta %>% select(sampleid, simpson_reciprocal), by = "sampleid") %>%
    mutate(pseudotiny = runif(n(), min = 0, max = 10^-10),
           changed_relab = relab + pseudotiny)

  keep <- prevalent_genera(g_relab, n_distinct(meta$sampleid))

  spearman_res <- g_relab %>%
    filter(genus %in% keep) %>%
    split(.$genus) %>%
    imap_dfr(function(.x, .y) {
      ct <- suppressWarnings(cor.test(.x$simpson_reciprocal, .x$changed_relab,
                                      method = "spearman", exact = TRUE))
      list(genus = .y, rho = ct$estimate, pval = ct$p.value)
    })

  spearman_res %>%
    mutate(padj = p.adjust(pval, method = "BH"),
           sig05 = if_else(padj < 0.05, "FDR < 0.05", "FDR >= 0.05"),
           Correlation = factor(if_else(rho >= 0, "higher_div", "lower_div"),
                                levels = c("lower_div", "higher_div"))) %>%
    arrange(rho, desc(genus))
}
