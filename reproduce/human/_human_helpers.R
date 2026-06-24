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
