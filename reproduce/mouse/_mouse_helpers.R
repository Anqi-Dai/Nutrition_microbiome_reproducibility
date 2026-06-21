# Shared helpers for the mouse figure panels.
#
# Every mouse panel reads a released input (the source-data workbook for CFU and
# behaviour panels, the long 16S relab table, or the RNA-seq count matrix) and
# renders one of a few recurring archetypes. The shared pieces, palettes, theme,
# log axis, statistical-comparison defaults and the save routine, live here so
# the panel scripts stay close to the originals (R01, R07, R08, R45, R38) while
# the styling is defined in exactly one place.

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggpubr)     # stat_compare_means: the in-figure significance brackets
  library(scales)     # log10 scientific-notation axis labels
  library(here)
})

# One data root, resolved from one environment variable. A public clone with no
# setup falls back to released_data; pointing NUTRITION_DATA at the full private
# folder runs the same code against the complete data with no path edits.
nutrition_data_root <- function() {
  Sys.getenv("NUTRITION_DATA", unset = here::here("released_data"))
}

mouse_workbook <- function() {
  file.path(nutrition_data_root(), "Dai_mouse_figure_raw_data.xlsx")
}

# All sequencing inputs are consolidated into the workbook: the RNA-seq gene count
# matrix as `raw_counts_matrix`, and the 16S relab table split across
# `mouse_16s_asv_relab` (Taxon x sample) and `mouse_16s_sample_meta` (per sample).

# Expensive fits (MaAsLin2, DESeq2, fgsea) cache here so plot reruns never refit.
intermediate_dir <- function() here::here("intermediate_data")
cache_path <- function(...) file.path(intermediate_dir(), ...)

# Sheet name is the panel-to-data map. Reading goes through one function so a
# missing workbook produces one clear, named error rather than a deep stack.
read_mouse_sheet <- function(sheet) {
  path <- mouse_workbook()
  if (!file.exists(path)) {
    stop("Mouse workbook not found at: ", path,
         "\nSet NUTRITION_DATA or place Dai_mouse_figure_raw_data.xlsx in released_data/.")
  }
  read_excel(path, sheet = sheet)
}

# The combined "Treatment + Day" label always ends in the day number, whatever
# separator the treatment name happens to use. Pulling the trailing integer is
# more robust than splitting on a separator that varies between sheets.
add_day <- function(df, from_col, to_col = "day") {
  df |>
    mutate("{to_col}" := as.integer(str_extract(.data[[from_col]], "\\d+$")))
}

# Manuscript palettes. The sucrose experiments use grey for vehicle and pink for
# sucrose, darker shades for the antibiotic arm; the smoothie experiment swaps
# pink for orange; the chow panels are antibiotic-only (two arms).
pal_sucrose4 <- c("gray76", "#ffbcdc", "gray32", "deeppink2")  # PBS veh, PBS suc, abx veh, abx suc
pal_smoothie4 <- c("gray76", "bisque2", "gray32", "darkorange3")
pal_chow2 <- c("gray32", "deeppink2")                          # abx veh, abx suc
# Delayed-sucrose experiment (R31): antibiotic-only vehicle dark grey, the sucrose
# timing arms pink, the PBS controls in light shades. Drives E8j (AUC) and E8k.
pal_delay6 <- c(
  "PBS + Water"              = "lightgray",
  "PBS + Sucrose"            = "lightpink",
  "Abx + Water"              = "gray32",
  "Abx + Sucrose (Standard)" = "deeppink2",
  "Abx + Sucrose (Delay +1)" = "deeppink2",
  "Abx + Sucrose (Delay +2)" = "deeppink2")
# Diet palette for the weight overlay (R31): regular water grey, sucrose pink.
pal_diet2 <- c(regular_water = "gray32", sucrose_water = "deeppink2")
# Volcano and lollipop direction colours, matching R01 / R38.
pal_volcano <- c(neg = "blue", not = "black", pos = "red")
pal_lollipop <- c(Downregulated = "royalblue", Upregulated = "#E64B35")

# Order the day-faceted x axis as group-major, day-minor (matches R01/R07): each
# treatment block sits together, days ascending within it.
xvar_levels <- function(grp_levels, days) {
  as.vector(t(outer(grp_levels, days, paste, sep = "__")))
}

# CFU axes are read on a log scale with 10^n ticks, as in every original CFU panel.
scale_log10_sci <- function() {
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))
}

# Several CFU panels encode the treatment group by colour while the x axis shows
# only the day, so the colour mapping is unreadable on its own. These render the
# group names as a right-side legend, prettifying the factor labels for display
# (PBS__vehicle -> PBS + vehicle) without touching the underlying data. Pair
# pretty_grp (passed to scale_colour_manual's labels) with legend_treatment (the
# theme override that re-enables the legend theme_mouse() switches off).
pretty_grp <- function(x) str_replace_all(x, "_+", " + ")

legend_treatment <- function() {
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.9, "lines"))
}

theme_mouse <- function(base_size = 11) {
  theme_light(base_size = base_size) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(colour = "black", face = "bold"),
      axis.text = element_text(colour = "black")
    )
}

# Panel-level outputs with manuscript-matching names, one file per panel, all
# landing in results/. The directory is created on demand so a fresh clone works.
save_panel <- function(plot, file, width = 3.4, height = 3.2) {
  out_dir <- here::here("results")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  path <- file.path(out_dir, file)
  ggsave(path, plot = plot, width = width, height = height, units = "in")
  message("wrote ", path)
  invisible(path)
}
