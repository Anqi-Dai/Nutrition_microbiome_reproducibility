# Fig 3 + Extended Data E6 clinical-outcome panels and Supplementary Tables 1-6
# (RESTRICTED). Ported from R09_clinical_outcome__code_for_Figure_3.Rmd.
#
# Consumes the cleaned, merged df_main (restricted_data/df_main_clinical_outcome.rds).
# That table is built upstream in the dev repo (Nutrition_microbiome/scripts/
# 62_build_clinical_df_main.R) and shipped here only in its cleaned form; this repo
# does not carry the build step.
# Produces:
#   F3a  forest: HR per abx-exposure day, above- vs below-median sugar density
#   F3b  forest: HR per abx-exposure day, diet-pattern Cluster 1 vs Cluster 2
#   F3   combined a/b
#   E6c  scatter: avg sugar intake (g) vs avg daily calories
#   E6d  scatter: avg sugar density (g/1000 kcal) vs avg daily calories
#   E6j  contour: log-HR over sugar density x abx exposure (spline Cox model)
#   Supplementary Tables 1-6 as one PDF (S1 full cohort, S2 microbiome
#        sub-cohort, S3/S5 characteristics by sugar / cluster, S4/S6 Cox models),
#        laid out like the accepted supplementary file: a guide page, a page of
#        legends, then one page per table with its title above and legend below
#
# Clinical outcomes are protected: reads the restricted df_main and skips cleanly
# when it is absent. (The E6i hospital-discharge cumulative-incidence panel is not
# here: it comes from R10, not R09.)

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({
  library(survival)
  library(gtsummary)
  library(patchwork)
  library(splines)
  library(grid)   # the supplementary tables are drawn as PDF pages, not a workbook
})

df_file <- "df_main_clinical_outcome.rds"
if (!has_restricted(df_file)) {
  message("Fig 3 / E6 skipped: restricted df_main not found (", restricted(df_file), ").")
  message("This cleaned table is built upstream in the dev repo; place it in restricted_data/.")
  quit(save = "no", status = 0)
}

df_main <- read_rds(restricted(df_file))
results_dir <- here::here("results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Cox models (landmarked at day 12; drop the one patient with OStime_30 == 0)
# ---------------------------------------------------------------------------
# NB: OStime_30 is a MISNOMER. It is defined upstream as `OStime - 12` (not - 30),
# so it already IS the day-12 landmark: filter(OStime_30 > 0) keeps the 172 patients
# who survived past day 12, matching the Fig. 3 legend. The name is legacy; there is
# no separate OStime_12 column. Same for tevent_30 (= tevent - 12).
fit_cluster <- coxph(
  Surv(OStime_30, OSevent) ~ modal_diet * day_exposed + intensity + source_and_gvhdppx,
  data = df_main |> filter(OStime_30 > 0))

fit_sugar_density_binarized <- coxph(
  Surv(OStime_30, OSevent) ~ SugarCal_cat_high * day_exposed + source_and_gvhdppx + intensity,
  data = df_main |> filter(OStime_30 > 0))

# ---------------------------------------------------------------------------
# Forests: per-group HR for abx exposure. The reference group's slope is
# `day_exposed`; the other group's slope is `day_exposed` + the interaction term,
# with SE from the covariance of the two.
# ---------------------------------------------------------------------------
group_forest <- function(fit, int_term, top_label, bottom_label, plot_title = NULL) {
  vc <- vcov(fit); b <- coef(fit)
  ref_est <- b["day_exposed"]; ref_se <- sqrt(vc["day_exposed", "day_exposed"])
  alt_est <- b["day_exposed"] + b[int_term]
  alt_se  <- sqrt(vc["day_exposed", "day_exposed"] + vc[int_term, int_term] +
                  2 * vc["day_exposed", int_term])
  tibble(group = c(top_label, bottom_label),
         estimate = c(ref_est, alt_est), se = c(ref_se, alt_se)) |>
    mutate(HR = exp(estimate),
           HR.low = exp(estimate - 1.96 * se), HR.high = exp(estimate + 1.96 * se),
           group = fct_inorder(group) |> fct_rev()) |>
    ggplot(aes(x = HR, y = group)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = HR.low, xmax = HR.high), height = 0, linewidth = 0.5) +
    geom_point(shape = 22, fill = "black", color = "black", size = 2) +
    scale_x_log10() +
    labs(x = "HR per day\nof antibiotic exposure", y = NULL, title = plot_title) +
    theme_classic(base_size = 12)
}

forest_sugar <- group_forest(
  fit_sugar_density_binarized, "SugarCal_cat_highBelow-median:day_exposed",
  "above-median\nsugar intake\ng/1000 kcal", "below-median\nsugar intake\ng/1000 kcal")

forest_cluster <- group_forest(
  fit_cluster, grep("modal_diet.*:day_exposed", names(coef(fit_cluster)), value = TRUE),
  "Cluster 1", "Cluster 2")

ggsave(file.path(results_dir, "F3a_sugar_density_HR_forest.pdf"), forest_sugar, width = 3, height = 2)
ggsave(file.path(results_dir, "F3b_cluster_HR_forest.pdf"), forest_cluster, width = 3, height = 2)

# report the headline numbers for verification
sugar_tab <- broom::tidy(fit_sugar_density_binarized, exponentiate = TRUE, conf.int = TRUE)
cluster_tab <- broom::tidy(fit_cluster, exponentiate = TRUE, conf.int = TRUE)
message("F3 check: sugar day_exposed HR = ",
        round(sugar_tab$estimate[sugar_tab$term == "day_exposed"], 2),
        "; cluster day_exposed HR = ",
        round(cluster_tab$estimate[cluster_tab$term == "day_exposed"], 2))

# ---------------------------------------------------------------------------
# E6c / E6d scatters (Spearman), sugar amount and density vs total calories
# ---------------------------------------------------------------------------
e6c <- ggplot(df_main, aes(avg_total_caloric_intake / 1000, avg_sugar_daily_intake_gram)) +
  geom_point(alpha = 0.5, color = "grey20") +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) +
  ggpubr::stat_cor(method = "spearman", label.x.npc = 0.02, label.y.npc = 0.98,
                   p.accuracy = 0.001, cor.coef.name = "R") +
  labs(x = expression("Average Daily Calories (" * 10^3 ~ "kcal)"),
       y = "Average Sugar Intake (g)") +
  theme_minimal(base_size = 12)

e6d <- ggplot(df_main, aes(avg_total_caloric_intake / 1000, avg_sugar_density_per_1000kcal)) +
  geom_point(alpha = 0.5, color = "grey20") +
  geom_smooth(method = "lm", se = FALSE, color = "darkred", formula = y ~ x) +
  ggpubr::stat_cor(method = "spearman", label.x.npc = 0.42, label.y.npc = 0.98,
                   p.accuracy = 0.001, cor.coef.name = "R") +
  labs(x = expression("Average Daily Calories (" * 10^3 ~ "kcal)"),
       y = "Average Sugar Proportion\n(g/1000 kcal)") +
  theme_minimal(base_size = 12)

ggsave(file.path(results_dir, "E6c_sugar_vs_calories.pdf"), e6c, width = 4, height = 3.2)
ggsave(file.path(results_dir, "E6d_sugardensity_vs_calories.pdf"), e6d, width = 4, height = 3.2)

# ---------------------------------------------------------------------------
# E6j contour: log-HR landscape over sugar density x abx exposure.
# Natural-spline sugar term x day_exposed, TCD+TCD / nonablative as reference.
# ---------------------------------------------------------------------------
df_use <- df_main |> filter(OStime_30 > 0)
knots <- quantile(df_use$avg_sugar_density_per_1000kcal, probs = c(0.1, 0.9))

fit_contour <- coxph(
  Surv(OStime_30, OSevent) ~
    ns(avg_sugar_density_per_1000kcal, knots = knots) * day_exposed +
    source_and_gvhdppx + intensity,
  data = df_use)

grid_df <- expand.grid(
  day_exposed = seq(0, 20, length.out = 120),
  avg_sugar_density_per_1000kcal = seq(23, 150, length.out = 120),
  source_and_gvhdppx = "TCD+TCD", intensity = "nonablative")

pred <- predict(fit_contour, newdata = grid_df, type = "lp", se.fit = TRUE)
grid_df <- grid_df |>
  mutate(lp = pred$fit, se = pred$se.fit, hr = exp(lp),
         lower_hr = exp(lp - 1.96 * se), upper_hr = exp(lp + 1.96 * se))

# diverging HR bins exactly as the published legend (14 bins, 0.08 .. 18.91). Clamp
# HR into that range so the extreme regions still land in the end bins rather than
# leaving white gaps. The bin straddling HR = 1 (0.78, 1.28] is the neutral (white)
# tile; the 6 bins below are blue, the 7 above are red, deepening to the extremes.
hr_breaks <- c(0.08, 0.11, 0.17, 0.25, 0.36, 0.53, 0.78, 1.28, 1.89, 2.77,
               4.06, 5.97, 8.77, 12.88, 18.91)
hr_fill <- c(colorRampPalette(c("#08306B", "#C6DBEF"))(6), "#F7F7F7",
             colorRampPalette(c("#FCBBA1", "#67000D"))(7))
grid_df <- grid_df |> mutate(hr = pmin(pmax(hr, 0.081), 18.9))
sugar_median <- median(df_use$avg_sugar_density_per_1000kcal, na.rm = TRUE)

main_contour <- ggplot(grid_df, aes(day_exposed, avg_sugar_density_per_1000kcal)) +
  geom_contour_filled(aes(z = hr), breaks = hr_breaks) +
  # faint grey dashed reference contour where HR = 1 (no significance mesh)
  geom_contour(aes(z = hr), breaks = 1, color = "grey70", linewidth = 1.1, linetype = 2) +
  geom_hline(yintercept = sugar_median, colour = "darkslateblue", linewidth = 0.8) +
  scale_fill_manual(values = hr_fill, drop = FALSE, name = "hazard ratio",
                    guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(23, 150)) +
  labs(x = "Days of broad-spectrum antibiotic exposure",
       y = "Average grams of sugar intake\nper 1000 Kcal between day -7 and day 12") +
  theme_classic(base_size = 11) +
  theme(legend.key.size = unit(0.35, "cm"))

# marginal patient histograms (top = abx exposure, right = sugar density)
top_hist <- ggplot(df_use, aes(day_exposed)) +
  geom_histogram(bins = 30, fill = "grey75", colour = "white", linewidth = 0.2) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 20)) +
  labs(y = "patients") +
  theme_classic(base_size = 9) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

right_hist <- ggplot(df_use, aes(avg_sugar_density_per_1000kcal)) +
  geom_histogram(bins = 30, fill = "grey75", colour = "white", linewidth = 0.2) +
  scale_x_continuous(expand = c(0, 0), limits = c(23, 150)) +
  labs(y = "patients") +
  coord_flip() +
  theme_classic(base_size = 9) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

e6j <- (top_hist + patchwork::plot_spacer() +
        main_contour + right_hist +
        plot_layout(widths = c(4, 1), heights = c(1, 4)))
ggsave(file.path(results_dir, "E6j_sugar_abx_HR_contour.pdf"), e6j, width = 8, height = 6)

# ---------------------------------------------------------------------------
# Supplementary Tables 1-6
# ---------------------------------------------------------------------------
# Shared clinical recode (labels, frequency ordering) used by the summary tables.
recode_char <- function(d) {
  d |>
    mutate(Disease = case_when(
      disease.simple == "NHL" ~ "Non-Hodgkin's lymphoma",
      disease.simple == "MDS/MPN" ~ "MDS/MPN",
      disease.simple == "AML" ~ "Acute myeloid leukemia",
      disease.simple == "ALL" ~ "Acute lymphoid leukemia",
      disease.simple == "CLL" ~ "Chronic lymphocytic leukemia",
      disease.simple %in% c("CML", "Hodgkins", "AA", "other") ~ "Other",
      disease.simple == "Myeloma" ~ "Myeloma"),
      source = case_when(
        source == "Unmodified" ~ "Unmodified bone marrow or PBSC",
        source == "Cord" ~ "Cord blood",
        source == "TCD" ~ "T-cell depleted PBSC"),
      sex = case_when(sex == "F" ~ "Female", sex == "M" ~ "Male"),
      intensity = case_when(
        intensity == "nonablative" ~ "Nonmyeloablative",
        intensity == "ablative" ~ "Ablative",
        intensity == "reduced" ~ "Reduced intensity"),
      gvhd_ppx = case_when(
        gvhd_ppx == "CNI-based" ~ "CNI-based",
        gvhd_ppx == "PTCy-based" ~ "PTCy-based",
        gvhd_ppx == "TCD" ~ "T-cell depleted PBSC")) |>
    mutate(intensity = fct_reorder(intensity, intensity, .fun = length, .desc = TRUE),
           source    = fct_reorder(source, source, .fun = length, .desc = TRUE),
           Disease   = fct_reorder(Disease, Disease, .fun = length, .desc = TRUE),
           sex       = fct_reorder(sex, sex, .fun = length, .desc = TRUE)) |>
    rename(`Graft type` = source,
           `Intensity of conditioning regimen` = intensity,
           Sex = sex, Age = age, `GvHD prophylaxis` = gvhd_ppx,
           `Days exposed to broad-spectrum antibiotics` = day_exposed,
           `Length of stay at hospital` = leng_of_stay)
}

# S1 / S2: overall characteristics, mean (SD). S2 restricted to the microbiome
# sub-cohort (patients with a stool sample in the released META).
overall_summary <- function(d) {
  recode_char(d) |>
    select(Age, Sex, Disease, `GvHD prophylaxis`,
           `Days exposed to broad-spectrum antibiotics`, `Length of stay at hospital`,
           `Graft type`, `Intensity of conditioning regimen`) |>
    tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})",
                                 all_categorical() ~ "{n} ({p}%)")) |>
    bold_labels()
}

micro_pids <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) |>
  distinct(pid) |> pull(pid)

s1 <- overall_summary(df_main)
s2 <- overall_summary(df_main |> filter(pid %in% micro_pids))

# S3 / S5: characteristics by group, median (IQR).
grouped_summary <- function(d, by_var) {
  recode_char(d) |>
    select(Age, Sex, Disease, `Graft type`, `Intensity of conditioning regimen`,
           `GvHD prophylaxis`, `Days exposed to broad-spectrum antibiotics`,
           all_of(by_var)) |>
    tbl_summary(by = all_of(by_var),
                statistic = list(all_continuous() ~ "{median} ({IQR})",
                                 all_categorical() ~ "{n} ({p}%)")) |>
    bold_labels()
}
s3 <- grouped_summary(df_main, "SugarCal_cat_high")
s5 <- grouped_summary(df_main, "modal_diet")

# S4 / S6: multivariable Cox models.
s4 <- tbl_regression(
  fit_sugar_density_binarized,
  include = c(SugarCal_cat_high, day_exposed, "SugarCal_cat_high:day_exposed"),
  exponentiate = TRUE,
  label = list(SugarCal_cat_high ~ "Averaged Sugar Intake",
               day_exposed ~ "Abx exposure (days)")) |>
  bold_p(t = 0.05) |> bold_labels()

s6 <- tbl_regression(
  fit_cluster,
  include = c(modal_diet, day_exposed, "modal_diet:day_exposed"),
  exponentiate = TRUE,
  label = list(modal_diet ~ "Diet-pattern cluster",
               day_exposed ~ "Abx exposure (days)")) |>
  bold_p(t = 0.05) |> bold_labels()

# ---------------------------------------------------------------------------
# PDF renderer for the supplementary tables
# ---------------------------------------------------------------------------
# The accepted supplementary file is a PDF rather than a workbook, so the tables
# are drawn straight onto PDF pages here instead of being written to sheets and
# pasted into a document by hand. The layout follows that file: US Letter
# landscape, a shaded header band closed by a single black rule and no other
# rules, bold variable labels, the title above each table and the full legend
# below it. All geometry is in points measured from the TOP-LEFT of the page, the
# way it was measured off the accepted file, and grid does the flip to its own
# bottom-left origin in one place (`put_text`).

PAGE_W <- 792   # 11 in
PAGE_H <- 612   # 8.5 in

# Calibri sets the accepted file but is a Microsoft font that is not installed
# everywhere; fall through its metric-compatible clone to a generic sans so the
# script runs on any machine. Only the shapes change, never the layout.
supp_font <- {
  installed <- if (requireNamespace("systemfonts", quietly = TRUE)) {
    unique(systemfonts::system_fonts()$family)
  } else character()
  found <- intersect(c("Calibri", "Carlito"), installed)
  if (length(found) > 0) found[1] else "Helvetica"
}

# Width of a string on the open device, in points. Everything that wraps or gets
# a column width is measured rather than guessed, so a longer label or an extra
# group column widens the table instead of overprinting the next one.
text_w <- function(s, size, face = "plain") {
  convertWidth(
    grobWidth(textGrob(s, gp = gpar(fontfamily = supp_font, fontface = face,
                                    fontsize = size))),
    "pt", valueOnly = TRUE)
}

# One string with its top-left corner at (x, y) from the top-left of the page.
put_text <- function(s, x, y, size, face = "plain") {
  grid.text(s, x = unit(x, "pt"), y = unit(PAGE_H - y, "pt"),
            just = c("left", "top"),
            gp = gpar(fontfamily = supp_font, fontface = face, fontsize = size))
}

# A paragraph whose leading run (the table number) is bold and whose remainder is
# regular, wrapped greedily at `width`. Words are placed one at a time because
# the bold run ends mid-line and the two faces have different widths. Returns the
# top of the line after the paragraph, so blocks stack without hardcoded heights.
put_paragraph <- function(bold_run, plain_run, x, y, width, size, leading) {
  split_words <- function(s) if (nchar(trimws(s)) == 0) character() else strsplit(trimws(s), " +")[[1]]
  words <- c(split_words(bold_run), split_words(plain_run))
  faces <- c(rep("bold", length(split_words(bold_run))),
             rep("plain", length(split_words(plain_run))))

  space_w <- text_w(" ", size)
  widths  <- vapply(seq_along(words), \(i) text_w(words[i], size, faces[i]), numeric(1))

  cur_x <- x
  line_y <- y
  for (i in seq_along(words)) {
    # break before a word that would run past the right edge, unless it is the
    # first word on the line (nothing to gain by breaking again)
    if (cur_x > x && cur_x + widths[i] > x + width) {
      cur_x <- x
      line_y <- line_y + leading
    }
    put_text(words[i], cur_x, line_y, size, faces[i])
    cur_x <- cur_x + widths[i] + space_w
  }
  line_y + leading
}

# gtsummary marks bold cells with markdown underscores (__cell__) and bold column
# headers with asterisks (**header**). Split a table into its plain text and a
# parallel logical matrix of which cells were bold, then drop the markers.
gt_cells <- function(t) {
  raw <- as_tibble(t)
  txt <- lapply(raw, \(col) ifelse(is.na(col), "", as.character(col)))

  bold <- vapply(txt, \(col) str_detect(col, "^__.*__$"), logical(length(txt[[1]])))
  bold <- matrix(bold, nrow = length(txt[[1]]), ncol = length(txt))

  txt <- lapply(txt, \(col) str_remove_all(col, "\\*\\*|__"))
  # group headers arrive as "**Above-median**  \nN = 86"; dropping the newline
  # (not the spaces around it) gives the single-line header the accepted file uses
  headers <- names(raw) |> str_remove_all("\\*\\*|__") |> str_replace_all("\n", "")
  colnames(bold) <- headers

  # In the accepted file a Cox row whose p-value cleared bold_p() is bold right
  # across, not only in the p-value cell. Carry that emphasis over.
  if ("p-value" %in% headers) bold[bold[, "p-value"], ] <- TRUE

  txt <- as.data.frame(txt, check.names = FALSE, stringsAsFactors = FALSE)
  names(txt) <- headers
  list(text = txt, bold = bold)
}

# The same structure for a table that never came from gtsummary: header bold,
# no bold body cells. Supplementary Table 7 is the only one of these.
plain_cells <- function(df) {
  txt <- as.data.frame(lapply(df, as.character), check.names = FALSE,
                       stringsAsFactors = FALSE)
  names(txt) <- names(df)
  list(text = txt,
       bold = matrix(FALSE, nrow = nrow(txt), ncol = ncol(txt),
                     dimnames = list(NULL, names(txt))))
}

# The abbreviation footnote under the Cox tables, read off the gtsummary object
# rather than retyped, restricted to the columns the table actually shows (the
# standard-error abbreviation is carried but never printed).
gt_abbrev <- function(t) {
  ab <- t$table_styling$abbreviation
  shown <- t$table_styling$header |> filter(!hide) |> pull(column)
  ab <- ab |> filter(column %in% shown) |> pull(abbreviation) |> unique() |> sort()
  if (length(ab) == 0) NULL else paste0("Abbreviations: ", paste(ab, collapse = ", "))
}

# Draw a table with its top-left cell corner at (x, y). Returns the y below the
# last row. Column widths come from the widest cell in each column, always
# measured in bold so a bolded cell never overflows its column.
put_table <- function(cells, x, y, size = 11, row_h = 14, head_h = 15,
                      body_gap = 3, pad = 28, min_width = 560) {
  txt <- cells$text
  bold <- cells$bold
  n_col <- ncol(txt)

  head_w <- vapply(names(txt), text_w, numeric(1), size = size, face = "bold")
  body_w <- vapply(txt, \(col) max(c(0, vapply(col, text_w, numeric(1),
                                               size = size, face = "bold"))),
                   numeric(1))
  natural <- pmax(head_w, body_w) + pad

  # value columns share one width so the numbers sit in even columns; the label
  # column absorbs any slack needed to reach `min_width`, which is what gives the
  # accepted file its wide gap between the labels and the first value column
  value_w <- if (n_col > 1) max(natural[-1]) else 0
  label_w <- natural[1]
  total_w <- label_w + value_w * (n_col - 1)
  if (total_w < min_width) {
    label_w <- label_w + (min_width - total_w)
    total_w <- min_width
  }
  col_x <- x + cumsum(c(0, label_w, rep(value_w, max(n_col - 2, 0))))

  # shaded header band, inset 3 pt to the left of the text as in the accepted
  # file, closed by the single black rule that is the table's only rule
  grid.rect(x = unit(x - 3, "pt"), y = unit(PAGE_H - y, "pt"),
            width = unit(total_w, "pt"), height = unit(head_h, "pt"),
            just = c("left", "top"), gp = gpar(fill = "#D9D9D9", col = NA))
  grid.rect(x = unit(x - 3, "pt"), y = unit(PAGE_H - (y + head_h - 1), "pt"),
            width = unit(total_w, "pt"), height = unit(1, "pt"),
            just = c("left", "top"), gp = gpar(fill = "black", col = NA))

  for (j in seq_len(n_col)) put_text(names(txt)[j], col_x[j], y + 1, size, "bold")

  # body_gap keeps the first row clear of the rule that closes the header band:
  # text is placed by the top of its ascent, so without it the first row would
  # sit right on the rule instead of below it
  for (i in seq_len(nrow(txt))) {
    row_y <- y + head_h + body_gap + (i - 1) * row_h
    for (j in seq_len(n_col)) {
      if (nchar(txt[[j]][i]) == 0) next
      put_text(txt[[j]][i], col_x[j], row_y, size,
               if (bold[i, j]) "bold" else "plain")
    }
  }
  list(y = y + head_h + body_gap + nrow(txt) * row_h, width = total_w)
}

# ---------------------------------------------------------------------------
# Supplementary Table 7, folded in from the mouse side
# ---------------------------------------------------------------------------
# S7 (the composition of the two mouse diets behind E8l) is NOT built here.
# `51_extdata_8.R` builds it next to the panel it belongs to, reading the control
# column off the LabDiet 5053 vendor data sheet rather than typing it in, and
# caches the table together with its title and legend to
# intermediate_data/Supplementary_Table_7_diet_composition.rds. This script only
# reads that cache, so the diet numbers and their wording keep exactly one owner:
# edit them in 51, not here.
#
# Run 51 before this script to get the full seven-table document. Without the
# cache the PDF is still written, with tables 1-6 and a guide page listing six,
# and the message at the end says which one came out.
supp7_cache <- cache_path("Supplementary_Table_7_diet_composition.rds")

supp7_entry <- if (file.exists(supp7_cache)) {
  s7 <- readRDS(supp7_cache)
  s7_number <- "Supplementary Table 7."
  # The cached object carries the number inside its title, and inside its legend
  # in markdown bold, because 51 renders both as a single string. Split it back
  # off so S7 is laid out with the same bold-number prefix as the other six.
  list(cells = plain_cells(s7$table |>
                             rename(Characteristic = characteristic,
                                    `Control diet` = control_diet,
                                    `No fiber diet` = no_fiber_diet)),
       abbrev = NULL,
       number = s7_number,
       title  = str_squish(str_remove(s7$title, fixed(s7_number))),
       legend = str_squish(str_remove(str_remove_all(s7$legend, "\\*\\*"),
                                      fixed(s7_number))))
} else {
  NULL
}

# ---------------------------------------------------------------------------
# What goes on each page
# ---------------------------------------------------------------------------
# Each table carries its drawn cells plus three pieces of text, kept apart
# because they are not the same words. `number` prefixes both the guide line and
# the legend and is the only part set in bold; `title` is the short title used on
# the guide page and as the heading above the table; `legend` is the full caption
# printed on the legends page and again under the table. Wording for 1-6 is the
# accepted supplementary file's, verbatim, so a regenerated PDF can be dropped in
# without re-editing; 7's comes from 51's cache.
supp_tables <- list(
  list(cells = gt_cells(s1), abbrev = gt_abbrev(s1),
       number = "Supplementary Table 1.",
       title = "Patient characteristics of the full nutrition cohort.",
       legend = paste(
         "Patient characteristics of the full nutrition cohort. Data are presented as mean",
         "(standard deviation) for continuous variables and as count (percentage) for categorical",
         "variables. PBSC, peripheral blood stem cells. MDS/MPN, myelodysplastic",
         "syndromes/myeloproliferative neoplasms. GVHD, graft-versus-host disease. CNI,",
         "calcineurin inhibitor. PTCy, post-transplant cyclophosphamide.")),
  list(cells = gt_cells(s2), abbrev = gt_abbrev(s2),
       number = "Supplementary Table 2.",
       title = "Patient characteristics of the microbiome-nutrition sub-cohort.",
       legend = paste(
         "Patient characteristics of the microbiome-nutrition sub-cohort. Data are presented as mean",
         "(standard deviation) for continuous variables and as count (percentage) for categorical",
         "variables. PBSC, peripheral blood stem cells. MDS/MPN, myelodysplastic",
         "syndromes/myeloproliferative neoplasms. GVHD, graft-versus-host disease. CNI,",
         "calcineurin inhibitor. PTCy, post-transplant cyclophosphamide.")),
  list(cells = gt_cells(s3), abbrev = gt_abbrev(s3),
       number = "Supplementary Table 3.",
       title = paste("Patient characteristics by dichotomized averaged sugar intake,",
                     "normalized by total caloric intake."),
       legend = paste(
         "Patient characteristics in the cohort dichotomized at the median sugar intake",
         "(above-median vs below-median), normalized by total caloric intake and averaged across",
         "HCT days -7 to 12. Data are presented as median (interquartile range, IQR) for continuous",
         "variables and as count (percentage) for categorical variables. PBSC, peripheral blood stem",
         "cells. MDS/MPN, myelodysplastic syndromes / myeloproliferative neoplasms. GVHD,",
         "graft-versus-host disease. CNI, calcineurin inhibitor. PTCy, post-transplant",
         "cyclophosphamide.")),
  list(cells = gt_cells(s4), abbrev = gt_abbrev(s4),
       number = "Supplementary Table 4.",
       title = paste("Multivariable Cox proportional hazards model for overall survival (OS)",
                     "in dichotomized sugar-intake groups."),
       legend = paste(
         "Multivariable Cox proportional hazards model for overall survival (OS) in dichotomized",
         "sugar-intake groups. The table displays the adjusted hazard ratios (HR) and 95% confidence",
         "intervals (CI) for the association between various predictors and OS, landmarked at day 12",
         "after HCT. The model includes an interaction term between sugar-intake group and days of",
         "broad-spectrum antibiotic exposure (between day -7 and day 12), adjusting for conditioning",
         "regimen intensity, graft type, and GVHD prophylaxis. The row labeled",
         "“Abx exposure(days)” corresponds to the HR in the above-median sugar-intake group.")),
  list(cells = gt_cells(s5), abbrev = gt_abbrev(s5),
       number = "Supplementary Table 5.",
       title = "Patient characteristics by diet-pattern cluster.",
       legend = paste(
         "Patient characteristics by diet-pattern cluster. The clusters were defined by",
         "latent-trajectory class analysis without taking clinical variables into account. Data are",
         "presented as median (interquartile range, IQR) for continuous variables and as count",
         "(percentage) for categorical variables. PBSC, peripheral blood stem cells. MDS/MPN,",
         "myelodysplastic syndromes / myeloproliferative neoplasms. GVHD, graft-versus-host disease.",
         "CNI, calcineurin inhibitor. PTCy, post-transplant cyclophosphamide.")),
  list(cells = gt_cells(s6), abbrev = gt_abbrev(s6),
       number = "Supplementary Table 6.",
       title = paste("Multivariable Cox proportional hazards model for overall survival (OS)",
                     "in dietary-pattern clusters."),
       legend = paste(
         "Multivariable Cox proportional hazards model for overall survival (OS) in dietary-pattern",
         "clusters. The table displays the adjusted hazard ratios (HR) and 95% confidence intervals",
         "(CI) for the association between various predictors and OS, landmarked at day 12 after HCT.",
         "The model includes an interaction term between diet-pattern cluster and days of",
         "broad-spectrum antibiotic exposure (between day -7 and day 12), adjusting for conditioning",
         "regimen intensity, graft type, and GVHD prophylaxis. The row labeled",
         "“Abx exposure(days)” corresponds to the HR in the cluster 1."))
)
if (!is.null(supp7_entry)) supp_tables <- c(supp_tables, list(supp7_entry))

# ---------------------------------------------------------------------------
# Draw the pages
# ---------------------------------------------------------------------------
supp_pdf <- file.path(results_dir, if (is.null(supp7_entry)) {
  "Supplementary_Tables_1_6.pdf"
} else {
  "Supplementary_Tables_1_7.pdf"
})
# cairo_pdf rather than the base pdf device: the base device draws ASCII hyphens
# ("sub-cohort", "graft-versus-host") with the font's long minus glyph, and it
# does not carry the curly quotes in the Cox legends.
grDevices::cairo_pdf(supp_pdf, width = PAGE_W / 72, height = PAGE_H / 72,
                     onefile = TRUE)

# Page 1, the guide: every table listed by number and short title.
grid.newpage()
put_text("Supplemental Guide", 21, 36, 14, "bold")
guide_y <- 69
for (t in supp_tables) {
  guide_y <- put_paragraph(t$number, t$title, x = 21, y = guide_y,
                           width = 750, size = 11, leading = 17)
}

# Page 2, the full legends, one paragraph per table separated by a blank line.
grid.newpage()
put_text("Supplemental Table Legends", 36, 36, 11, "bold")
legend_y <- 36 + 2 * 13.45
for (t in supp_tables) {
  legend_y <- put_paragraph(t$number, t$legend, x = 36, y = legend_y,
                            width = 722, size = 11, leading = 13.45)
  legend_y <- legend_y + 13.45
}

# One page per table: heading, table, abbreviation footnote if the table has one,
# then the full legend. The legend is wrapped to the table's own width so the
# block reads as one unit however wide the table turned out.
for (t in supp_tables) {
  grid.newpage()
  heading_y <- put_paragraph(paste(t$number, t$title), "", x = 39, y = 36,
                             width = 717, size = 14, leading = 18)

  drawn <- put_table(t$cells, x = 39, y = heading_y + 13)

  foot_y <- drawn$y
  if (!is.null(t$abbrev)) {
    put_text(t$abbrev, 39, foot_y + 1, 10, "italic")
    foot_y <- foot_y + 14
  }
  put_paragraph(t$number, t$legend, x = 39, y = foot_y + 13,
                width = max(drawn$width, 520), size = 11, leading = 13.45)
}
invisible(grDevices::dev.off())

message("Wrote Fig 3 / E6 panels and ", basename(supp_pdf), " to results/ (",
        length(supp_tables), " tables, font: ", supp_font, ").")
if (is.null(supp7_entry)) {
  message("Supplementary Table 7 was left out: no ", basename(supp7_cache),
          " in intermediate_data/. Run reproduce/51_extdata_8.R first to include it.")
}
