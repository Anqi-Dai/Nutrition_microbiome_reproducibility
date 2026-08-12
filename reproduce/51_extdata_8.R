# Extended Data Figure E8, non-sequencing panels.
#
# E8 a-f come from 16S and are produced by 53_extdata_8_16s.R. This file covers
# the panels that read a self-contained source-data sheet and need no count table:
#   E8g  daily chow consumption per day (biapenem + sucrose vs + vehicle)
#   E8h  trapezoidal AUC of chow consumption
#   E8j  AUC of Enterococcus under delayed sugar, both experiments pooled
#   E8k  delayed-sucrose CFU/mg time course, one facet per treatment
#   E8l  Enterococcus CFU/mg on fibre-free chow, already in log10 units
#   E8m  body-weight % change from baseline, per-mouse overlay with median + IQR
#
# Style and statistical comparisons follow the originals R08 (E8g/h) and R31
# (E8j AUC, E8k longitudinal, E8m weight overlay) and R07 (E8l).

source(here::here("reproduce", "mouse", "_mouse_helpers.R"))

# Shared with E8h: the chow AUC sheet still labels the arms without spaces.
chow_levels <- c("biapenem\n+\nvehicle", "biapenem\n+\nsucrose")
chow_comp <- list(c("biapenem\n+\nvehicle", "biapenem\n+\nsucrose"))

# Delayed-sucrose design, shared by E8j (AUC) and E8k (time course). The plotting
# order puts the antibiotic-only arm first; E8k facets in reading order (its
# reverse), E8j stacks them bottom-to-top after the coordinate flip.
delay_order <- c("PBS + Sucrose", "PBS + Water", "Abx + Sucrose (Delay +2)",
                 "Abx + Sucrose (Delay +1)", "Abx + Sucrose (Standard)", "Abx + Water")

# E8g: daily chow consumed, per day (R08) -------------------------------------
# The updated sheet restores the per-day label (Day3/Day6/Day9), so this returns
# to R08's per-day facet rather than the earlier pooled boxplot. This sheet labels
# the arms with spaces around the separator, so the levels are matched from the
# data instead of a hard-coded string.
chow_day <- read_mouse_sheet("Sup_Figure_8g_chow_consumed_day") |>
  rename(chow = chow_per_day, grp = treatment) |>
  mutate(day = factor(day, levels = c("Day3", "Day6", "Day9")))

chow_g_levels <- c(grep("vehicle", unique(chow_day$grp), value = TRUE),
                   grep("sucrose", unique(chow_day$grp), value = TRUE))
chow_day <- chow_day |> mutate(grp = factor(grp, levels = chow_g_levels))

p_8g <- chow_day |>
  ggplot(aes(x = grp, y = chow, colour = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, shape = 16, size = 1) +
  facet_wrap(~day) +
  scale_colour_manual(values = pal_chow2) +
  scale_x_discrete(labels = c("vehicle", "sucrose")) +
  stat_compare_means(comparisons = list(chow_g_levels), label = "p.signif",
                     method = "wilcox.test", tip.length = 0.04) +
  labs(x = NULL, y = "Chow consumption (g)\nper mouse per day") +
  theme_mouse()

save_panel(p_8g, "E8g_chow_per_day.pdf", width = 4.6, height = 3.4)

# E8h: AUC of chow consumed (R08) ---------------------------------------------
chow_auc <- read_mouse_sheet("Sup_Figure_8h_chow_consumed_auc") |>
  rename(auc = `Trapezoidal AUC (grams of Chow Consumed)`, grp = `Treatment Group`) |>
  mutate(grp = factor(grp, levels = chow_levels))

p_8h <- chow_auc |>
  ggplot(aes(x = grp, y = auc, colour = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, shape = 16, size = 1) +
  scale_colour_manual(values = pal_chow2) +
  stat_compare_means(comparisons = chow_comp, label = "p.signif",
                     method = "wilcox.test", tip.length = 0.04) +
  labs(x = NULL, y = "Trapezoidal\nAUC") +
  theme_mouse()

save_panel(p_8h, "E8h_chow_auc.pdf", width = 2.6, height = 3.4)

# E8j: delayed-sugar AUC, both experiments pooled (R31) -----------------------
# The published panel pools the two experiments (no facet) into a horizontal
# layout, encodes the experiment by point shape, and reports Wilcoxon p-values
# comparing the antibiotic-only arm against each sucrose-timing arm.
delay_comp <- list(c("Abx + Water", "Abx + Sucrose (Standard)"),
                   c("Abx + Water", "Abx + Sucrose (Delay +1)"),
                   c("Abx + Water", "Abx + Sucrose (Delay +2)"))

delay_auc <- read_mouse_sheet("Sup_Figure_8j_delayed_sucrose_a") |>
  rename(auc = `AUC of log10 (CFU per mg stool) x days`, grp = `Treatment Group`,
         experiment = `Experiment Number`) |>
  mutate(grp = factor(grp, levels = delay_order),
         experiment = str_replace(experiment, "_", " "))

p_8j <- delay_auc |>
  ggplot(aes(x = grp, y = auc, colour = grp)) +
  geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.8) +
  geom_jitter(aes(shape = experiment), width = 0.18, size = 1.5, alpha = 0.9) +
  scale_shape_manual(values = c("Experiment 1" = 16, "Experiment 2" = 15)) +
  scale_colour_manual(values = pal_delay6) +
  stat_compare_means(comparisons = delay_comp, label = "p.format",
                     method = "wilcox.test", tip.length = 0.01) +
  coord_flip(ylim = c(0, max(delay_auc$auc, na.rm = TRUE) * 1.25)) +
  labs(x = NULL, y = "Area Under the Curve\n(log10 CFU/mg * days)") +
  theme_mouse() +
  theme(axis.title.y = element_blank())

save_panel(p_8j, "E8j_delayed_sugar_auc.pdf", width = 5.2, height = 3.6)

# E8k: delayed-sucrose CFU time course (R31) ----------------------------------
# Raw per-mouse long table, reproduced as R31's longitudinal plot: one facet per
# treatment, day boxplots with each mouse's trajectory as a thin grey line, the
# per-day median connected, and point shape encoding the experiment. The facet
# strips name each treatment group, so no colour legend is needed.
cfu_8k <- read_mouse_sheet("Sup_Figure_8k_delayed_sucrose_C") |>
  rename(grp = treatment, mouse = mouse_identifier, experiment = experiment_no) |>
  mutate(cfu = as.numeric(cfus_per_mg_stool),  # column ships as text with literal "NA"
         grp = factor(grp, levels = rev(delay_order)),
         day = factor(day),
         experiment = str_replace(experiment, "_", " ")) |>
  filter(!is.na(cfu))

p_8k <- cfu_8k |>
  ggplot(aes(x = day, y = cfu, colour = grp)) +
  facet_wrap(~grp, ncol = 6) +
  geom_boxplot(outlier.shape = NA, aes(group = interaction(day, grp))) +
  geom_line(aes(group = mouse), colour = "gray70", alpha = 0.4) +
  geom_point(aes(shape = experiment), size = 1.2) +
  scale_shape_manual(values = c("Experiment 1" = 16, "Experiment 2" = 15)) +
  stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 0.9) +
  scale_colour_manual(values = pal_delay6) +
  scale_log10_sci() +
  labs(x = "Day of experiment", y = "Enterococcal\nCFU/mg") +
  theme_mouse() +
  theme(strip.text = element_text(size = 7))

save_panel(p_8k, "E8k_cfu_over_days.pdf", width = 8.0, height = 3.0)

# E8l: fibre-free chow, log10 CFU (R07) ---------------------------------------
groups_dpbs <- c("DPBS__vehicle", "DPBS__sucrose", "biapenem__vehicle", "biapenem__sucrose")
days_8l <- c(1, 3, 6, 9)

cfu_8l <- read_mouse_sheet("Sup_Figure_8l_no_fiber_chow_CFU") |>
  rename(log_cfu = `log10 (Enterococcal CFU per mg stool)`, grp = `Treatment Group`) |>
  add_day("Treatment + Day") |>
  mutate(grp = factor(grp, levels = groups_dpbs),
         xvar = factor(str_glue("{grp}__{day}"), levels = xvar_levels(groups_dpbs, days_8l)))

# Already log10, so the y axis stays linear; R07 tests biapenem veh vs suc on day 9 with a t-test.
p_8l <- cfu_8l |>
  ggplot(aes(x = xvar, y = log_cfu, colour = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, shape = 16, size = 0.9) +
  scale_colour_manual(values = pal_sucrose4, name = "Treatment group", labels = pretty_grp) +
  scale_x_discrete(labels = rep(days_8l, length(groups_dpbs))) +
  stat_compare_means(
    comparisons = list(c("biapenem__vehicle__9", "biapenem__sucrose__9")),
    label = "p.signif", method = "t.test", tip.length = 0.015, step.increase = 0.1) +
  labs(x = "Day", y = expression(log[10] ~ "Enterococcal CFU/mg")) +
  theme_mouse() +
  legend_treatment()

# The sheet carries no mouse identifier, so the animal count comes from the design
# rather than a distinct count. Sampling is balanced (every mouse plated on all four
# days, no dropouts; checked against the R07 raw plating table by
# internal/accepted_in_principle/08_cfu_unit_correction/verify_cfu_raw_tables.R), so each arm's mouse count is its
# row count on any one day.
n_mice_8l <- cfu_8l |>
  count(grp, day) |>
  group_by(grp) |>
  summarise(mice = max(n), .groups = "drop") |>
  pull(mice) |>
  sum()

print(paste("There are ", n_mice_8l, "individual mice."))
print(paste("There are ", nrow(cfu_8l), "individual samples."))

save_panel(p_8l, "E8l_fiberfree_log10_cfu.pdf", width = 6.0, height = 3.2)

# Supplementary Table 7: composition of the two diets behind E8l --------------
# The E8l legend points at Supplementary Table 7 for the diets, so the table is
# rebuilt here next to the panel it belongs to. It is a reference table, not a
# published figure panel, so it lands in intermediate_data/ rather than results/
# (results/ holds one PDF per manuscript panel and nothing else).
#
# The two columns have different provenance, and only one of them has a public
# source, so they are built differently:
#
# Control diet = LabDiet PicoLab Rodent Diet 20, catalog #5053. Every value is a
# field on the vendor product data sheet, so the column is READ OFF THAT SHEET
# rather than typed in: the PDF is downloaded once, cached, and parsed below. The
# odd phrase "percent of ration" in the row labels is the data sheet's own
# convention ("Nutrients expressed as percent of ration"), which is why the
# published table uses it.
#
# No-fiber diet = a fiber-free AIN-93G, ordered through W.F. Fisher and Son.
# There is NO public data sheet for it. The catalog AIN-93G that W.F. Fisher
# distributes (Bio-Serv F3156 pellet / F3197 powder) still contains ~4.8%
# cellulose and is 3.74 kcal/g, so it cannot be the 0%-fiber, 4.10 kcal/g diet
# this column describes; the column has to be a custom cellulose-free build. That
# column is therefore carried through as published and sanity-checked further
# down against the AIN-93G recipe with its cellulose swapped for cornstarch.
suppressPackageStartupMessages(library(ggtext))  # markdown title and legend text

# Pull the LabDiet data sheet once and keep it, so reruns and offline runs work.
# Values have been stable across revisions (the 09/27/22 and 09/05/23 sheets carry
# the same five numbers), but the parse is anchored on the field labels rather
# than on line positions so a re-typeset revision does not silently shift a cell.
labdiet_5053_url <- "https://www.labdiet.com/getmedia/253eaa5f-bc52-47cf-8dad-fe5828dc64c7/5053.pdf?ext=.pdf"
labdiet_5053_pdf <- cache_path("labdiet_5053_datasheet.pdf")

if (!file.exists(labdiet_5053_pdf)) {
  if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)
  try(curl::curl_download(labdiet_5053_url, labdiet_5053_pdf, quiet = TRUE), silent = TRUE)
}

# Take the number that follows a field label on the data sheet. pdf_text keeps the
# sheet's three-column layout, so one text line can hold fields from all three
# columns: cut everything up to and including the label, then take the first
# number after it. The dot leaders ("Fiber (Crude), % . . . . 4.4") are separated
# by spaces, so they never look like part of the value.
# Labels carry regex metacharacters ("Fiber (Crude), %"), so the label is located
# as a fixed string and the line is cut by position rather than by a pattern.
labdiet_field <- function(lines, label, after = NULL) {
  if (!is.null(after)) lines <- lines[seq(grep(after, lines, fixed = TRUE)[1], length(lines))]
  hit <- lines[grep(label, lines, fixed = TRUE)[1]]
  at <- regexpr(label, hit, fixed = TRUE)
  rest <- substring(hit, at + attr(at, "match.length"))
  as.numeric(regmatches(rest, regexpr("[0-9]+(\\.[0-9]+)?", rest)))
}

# The published control column, kept as the fallback for an offline run and as the
# assertion target when the data sheet is available.
supp7_control_published <- c(4.11, 4.4, 24.495, 13.122, 62.382)

if (file.exists(labdiet_5053_pdf) && requireNamespace("pdftools", quietly = TRUE)) {
  sheet_lines <- unlist(strsplit(pdftools::pdf_text(labdiet_5053_pdf), "\n"))
  # "Protein, %" and "Fat (ether extract), %" each appear twice on the sheet, once
  # as a percent of ration and once as a percent of calories, so the three calorie
  # fields are read only from the "Calories provided by:" block onward.
  supp7_control <- c(
    labdiet_field(sheet_lines, "Gross Energy, kcal/gm"),
    labdiet_field(sheet_lines, "Fiber (Crude), %"),
    labdiet_field(sheet_lines, "Protein, %", after = "Calories provided by:"),
    labdiet_field(sheet_lines, "Fat (ether extract), %", after = "Calories provided by:"),
    labdiet_field(sheet_lines, "Carbohydrates, %", after = "Calories provided by:"))
  # The sheet revision date, recorded so a future mismatch can be traced to a
  # vendor reformulation rather than to a parsing failure.
  labdiet_5053_revision <- regmatches(
    sheet_lines, regexpr("[0-9]{2}/[0-9]{2}/[0-9]{2}", sheet_lines)) |> unlist() |> head(1)
  message("LabDiet 5053 data sheet, revision ", labdiet_5053_revision, ": ",
          paste(supp7_control, collapse = ", "))
  stopifnot(isTRUE(all.equal(supp7_control, supp7_control_published)))
} else {
  supp7_control <- supp7_control_published
  labdiet_5053_revision <- NA_character_
  message("LabDiet 5053 data sheet not available offline; using the published values.")
}

# Assemble the table. The control column is formatted back to the printed strings
# from the parsed numbers; the no-fiber column is the published text.
supp7 <- tibble::tibble(
  characteristic = c("Energy (Kcal/g)", "Fiber percent of ration",
                     "Calories provided by protein", "Calories provided by fat",
                     "Calories provided by carbohydrates"),
  control_diet = c(format(supp7_control[1], nsmall = 2),
                   paste0(supp7_control[2], "%"),
                   paste0(supp7_control[3:5], "%")),
  no_fiber_diet = c("4.10", "0", "17.4%", "15.4%", "67.2%"))

# Title and legend exactly as they read in the supplementary file, kept with the
# data so the saved object is self-describing.
supp7_title <- paste(
  "Supplementary Table 7. Nutritional information of standard chow and",
  "no-fiber diet used in mouse experiments.")

supp7_legend <- paste(
  "**Supplementary Table 7.** Nutritional information of standard chow and",
  "no-fiber diet used in mouse experiments. The table displays the nutritional",
  "composition of the two chow diets used in the mouse experiments. The standard",
  "mouse diet in the vivarium is a grain-based chow with 4.4% fiber content",
  "(LabDiet PicoLab Rodent Diet 20 catalog #5053). The no-fiber diet is a",
  "purified diet with 0% fiber content (W.F. Fisher and Son; catalog #AIN-93G).")

# Sanity check on the NO-FIBER column only. The control column is asserted against
# the vendor data sheet above and needs nothing further; this column has no data
# sheet to check against, so the closest available test is to rebuild it from the
# published AIN-93G formula (Reeves 1993) in g per kg, replace the 50 g of
# cellulose with cornstarch, and recompute energy and the calorie split with
# Atwater factors 4/9/4. Casein is 85-90% protein by weight; 87% is used here as
# the mid-point, and the L-cystine counts toward protein. The result is an
# approximation of an unknown vendor formulation, not a reproduction of it, so
# small differences are expected and it never feeds the table.
ain93g_fiber_free <- c(cornstarch = 397.486 + 50, casein = 200, maltodextrin = 132,
                       sucrose = 100, soybean_oil = 70, cellulose = 0,
                       mineral_mix = 35, vitamin_mix = 10, l_cystine = 3,
                       choline_bitartrate = 2.5)

supp7_check <- {
  g_protein <- ain93g_fiber_free[["casein"]] * 0.87 + ain93g_fiber_free[["l_cystine"]]
  g_fat     <- ain93g_fiber_free[["soybean_oil"]]
  g_carb    <- sum(ain93g_fiber_free[c("cornstarch", "maltodextrin", "sucrose")])
  kcal      <- c(protein = g_protein * 4, fat = g_fat * 9, carb = g_carb * 4)
  tibble::tibble(
    characteristic = supp7$characteristic,
    published      = c(4.10, 0, 17.4, 15.4, 67.2),
    recomputed     = c(sum(kcal) / 1000, 0, 100 * kcal / sum(kcal))) |>
    mutate(difference = round(recomputed - published, 2),
           recomputed = round(recomputed, 2))
}

message("Supplementary Table 7, no-fiber column vs an approximated cellulose-free AIN-93G ",
        "(approximation only, the control column is asserted against the vendor sheet):")
print(supp7_check)

# Same question from the other direction, and the more useful output of the two:
# invert the published calorie split back into an as-fed composition. The control
# column proves the table follows the LabDiet/TestDiet convention, where "calories
# provided by" is protein/fat/carbohydrate weight x 4/9/4 over the physiological
# fuel value, so dividing the published percentages by that energy gives the grams
# per 100 g the vendor must have formulated to. This is the description to quote
# when asking W.F. Fisher and Son or the manufacturer to identify the diet.
supp7_implied <- {
  kcal_per_100g <- 4.10 * 100
  tibble::tibble(
    component   = c("protein", "fat", "carbohydrate"),
    pct_calories = c(17.4, 15.4, 67.2),
    atwater      = c(4, 9, 4)) |>
    mutate(g_per_100g = round(kcal_per_100g * pct_calories / 100 / atwater, 1))
}

message("Implied as-fed composition of the no-fiber diet (g per 100 g), for a ",
        "vendor lookup. Compare with AIN-93G: casein 200 g/kg (~17.8% protein), ",
        "soybean oil 70 g/kg, carbohydrate 63% plus the 5% cellulose:")
print(supp7_implied)

# Render the table the way it is laid out in the supplementary file: a shaded
# header band, a rule under it, no other rules, the title above and the full
# legend below. Everything is placed on a blank 0-1 canvas so the geometry is
# explicit rather than inherited from a table package's defaults.
supp7_long <- supp7 |>
  rename(`Characteristic` = characteristic, `Control diet` = control_diet,
         `No fiber diet` = no_fiber_diet) |>
  mutate(row = row_number()) |>
  pivot_longer(-row, names_to = "column", values_to = "value")

# Column left edges and row centres, derived from the table size rather than
# hand-tuned, so adding a row does not need the coordinates re-eyeballed.
supp7_x <- c(`Characteristic` = 0.005, `Control diet` = 0.40, `No fiber diet` = 0.68)
supp7_row_h <- 0.075
supp7_head_y <- 0.80
supp7_body_y <- supp7_head_y - supp7_row_h * seq_len(nrow(supp7))

supp7_cells <- supp7_long |>
  mutate(x = supp7_x[column], y = supp7_body_y[row])

supp7_header <- tibble::tibble(column = names(supp7_x), x = supp7_x,
                               y = supp7_head_y)

p_supp7 <- ggplot() +
  # shaded header band plus the single black rule that closes it
  geom_rect(aes(xmin = 0, xmax = 1,
                ymin = supp7_head_y - supp7_row_h / 2,
                ymax = supp7_head_y + supp7_row_h / 2),
            fill = "#D9D9D9") +
  geom_segment(aes(x = 0, xend = 1,
                   y = supp7_head_y - supp7_row_h / 2,
                   yend = supp7_head_y - supp7_row_h / 2),
               linewidth = 0.6, colour = "black") +
  # header labels and body cells, all left aligned on their column edge
  geom_text(data = supp7_header, aes(x = x, y = y, label = column),
            hjust = 0, fontface = "bold", size = 3.5) +
  geom_text(data = supp7_cells, aes(x = x, y = y, label = value),
            hjust = 0, size = 3.5) +
  # title above the table and the wrapped legend below it
  geom_textbox(aes(x = 0, y = 0.99, label = supp7_title),
               hjust = 0, vjust = 1, halign = 0, width = unit(6.6, "in"),
               size = 3.9, fontface = "bold", box.colour = NA, fill = NA,
               box.padding = unit(rep(0, 4), "pt")) +
  geom_textbox(aes(x = 0, y = min(supp7_body_y) - 0.08, label = supp7_legend),
               hjust = 0, vjust = 1, halign = 0, width = unit(6.6, "in"),
               size = 3.2, box.colour = NA, fill = NA,
               box.padding = unit(rep(0, 4), "pt")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void()

# Save the object (data + title + legend + the consistency check) so the table
# can be re-rendered or written out elsewhere without re-reading this script,
# and a PDF alongside it for eyeballing against the supplementary file.
supp7_object <- list(table = supp7, title = supp7_title, legend = supp7_legend,
                     check = supp7_check, implied = supp7_implied, plot = p_supp7,
                     source = list(control_diet = labdiet_5053_url,
                                   control_diet_revision = labdiet_5053_revision,
                                   no_fiber_diet = "no public data sheet located"))

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)
saveRDS(supp7_object, cache_path("Supplementary_Table_7_diet_composition.rds"))
# cairo_pdf rather than the base pdf device: the base device draws the ASCII
# hyphen in "no-fiber" and "#AIN-93G" with the font's long minus glyph, which
# reads as an en dash.
ggsave(cache_path("Supplementary_Table_7_diet_composition.pdf"), plot = p_supp7,
       width = 7.0, height = 4.2, units = "in", device = grDevices::cairo_pdf)
message("wrote ", cache_path("Supplementary_Table_7_diet_composition.pdf"))

# E8m: median weight percent-change by treatment (4 groups) -------------------
# The revised sheet is already the analysis-ready per-mouse table (treatment_group
# labelled to match the palette, day, percent_change). Reduce to the per-group
# per-day median with the 25th-75th percentile as error bars, then draw the four
# dodged median lines, following the reference code for this panel.
df_plot_ready <- read_mouse_sheet("Sup_Figure_8m_all_mouse_weight_") |>
  mutate(percent_change = as.numeric(percent_change),
         day = as.integer(day),
         treatment_group = factor(treatment_group,
           levels = c("PBS + Water", "PBS + Sucrose", "Abx + Water", "Abx + Sucrose")))

summary_data <- df_plot_ready |>
  group_by(treatment_group, day) |>
  summarise(median_change = median(percent_change, na.rm = TRUE),
            q1 = quantile(percent_change, 0.25, na.rm = TRUE),
            q3 = quantile(percent_change, 0.75, na.rm = TRUE),
            .groups = "drop")

pal_weight4 <- c("PBS + Water" = "lightgray", "PBS + Sucrose" = "lightpink",
                 "Abx + Water" = "gray32", "Abx + Sucrose" = "deeppink2")
dodge_width <- position_dodge(width = 0.1)

p_8m <- summary_data |>
  ggplot(aes(x = day, y = median_change, group = treatment_group, colour = treatment_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_line(linewidth = 1.2, position = dodge_width) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.2, linewidth = 0.8, position = dodge_width) +
  geom_point(size = 3, position = dodge_width) +
  scale_colour_manual(values = pal_weight4) +
  scale_x_continuous(breaks = c(1, 3, 6)) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Day of experiment", y = "Weight Change from Baseline (%)",
       colour = "Treatment group") +
  theme_bw(base_size = 12) +
  theme(legend.title = element_text(face = "bold"), legend.position = "bottom")

save_panel(p_8m, "E8m_weight_change.pdf", width = 5.2, height = 4.6)
