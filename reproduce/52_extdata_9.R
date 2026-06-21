# Extended Data Figure E9, non-sequencing panels.
#
# E9 e, f, g come from RNA-seq and are produced by 54_extdata_9_rnaseq.R. This
# file covers the panels that read a self-contained source-data sheet:
#   E9a  prolonged 21-day CFU time course, antibiotic x sugar
#   E9b  commercial smoothie (mixed-sugar) exposure          (style from R45)
#   E9c  biapenem with different sugars (same source sheet as E9a)
#   E9d  germ-free monocolonization: per-mouse time course + AUC boxplot
#   E9h  CFU in bone-marrow-transplant recipients            (style from R02)
#
# E9a / E9c read one identical sheet; E9a shows the full antibiotic x sugar design
# and E9c the biapenem-only sugar comparison. The many-group time courses keep the
# palette and log axis but omit significance brackets (too many pairs to draw
# cleanly); the originals report those tests in the text.

source(here::here("reproduce", "mouse", "_mouse_helpers.R"))

sugar_order <- c("vehicle", "glucose", "fructose", "sucrose")
pal_sugar4 <- c(vehicle = "gray55", glucose = "#8da0cb", fructose = "#fc8d62", sucrose = "deeppink2")

# E9a and E9c now ship as separate sheets: E9a is the 21-day vehicle-vs-sucrose
# course (days 1-21), E9c the biapenem alternate-sugar comparison (days 1,3,6).
# Both encode antibiotic, sugar and day in one combined key column.
read_course <- function(sheet) {
  read_mouse_sheet(sheet) |>
    rename(key = `abx_treatment__diet_treatment__day`, cfu = cf_us_per_gram_stool) |>
    mutate(cfu = as.numeric(cfu)) |>   # the 9c sheet ships as text with literal "NA"
    filter(!is.na(cfu)) |>
    separate(key, into = c("antibiotic", "sugar", "day"), sep = "__", remove = FALSE) |>
    mutate(sugar = droplevels(factor(sugar, levels = sugar_order)),
           day = factor(as.integer(day)))
}

# E9a: 21-day time course (vehicle vs sucrose) --------------------------------
course_9a <- read_course("Sup_Figure_9a_21_day_exp")

p_9a <- course_9a |>
  ggplot(aes(x = day, y = cfu, colour = sugar)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
             size = 0.8, alpha = 0.5, shape = 16) +
  facet_wrap(~antibiotic) +
  scale_colour_manual(values = pal_sugar4) +
  scale_log10_sci() +
  labs(x = "Day", y = "Enterococcal\nCFU/gram") +
  theme_mouse() +
  theme(legend.position = "right", legend.title = element_blank())

save_panel(p_9a, "E9a_21day_timecourse.pdf", width = 6.6, height = 3.4)

# E9c: alternate sugars, full antibiotic x diet design ------------------------
# Both antibiotics (PBS in light shades, biapenem in saturated ones) crossed with
# the four diets, laid out group-major / day-minor with day on the x axis. The
# brackets compare the biapenem vehicle arm against each sugar at day 3 and day 6,
# the figure's key contrasts (vehicle bloom vs sugar-amplified bloom).
grp_9c <- c("PBS__vehicle", "PBS__sucrose", "PBS__glucose", "PBS__fructose",
            "biapenem__vehicle", "biapenem__sucrose", "biapenem__glucose", "biapenem__fructose")
pal_9c <- c("PBS__vehicle" = "gray76",  "PBS__sucrose" = "#ffbcdc",
            "PBS__glucose" = "#e2c290", "PBS__fructose" = "#aed9a0",
            "biapenem__vehicle" = "gray32",  "biapenem__sucrose" = "deeppink2",
            "biapenem__glucose" = "#8a6d2f", "biapenem__fructose" = "#3f8f3f")
days_9c <- c(1, 3, 6)

course_9c <- read_course("Sup_Figure_9c_alternate_sugars") |>
  mutate(antibiotic = recode(antibiotic, DPBS = "PBS"),
         day = as.integer(as.character(day)),
         grp = factor(str_glue("{antibiotic}__{sugar}"), levels = grp_9c),
         xvar = factor(str_glue("{antibiotic}__{sugar}__{day}"),
                       levels = xvar_levels(grp_9c, days_9c)))

p_9c <- course_9c |>
  ggplot(aes(x = xvar, y = cfu, colour = grp)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, shape = 16, size = 0.8) +
  scale_colour_manual(values = pal_9c, name = "Antibiotic + diet", labels = pretty_grp) +
  scale_log10_sci() +
  scale_x_discrete(labels = rep(days_9c, length(grp_9c))) +
  stat_compare_means(
    comparisons = list(c("biapenem__vehicle__3", "biapenem__sucrose__3"),
                       c("biapenem__vehicle__6", "biapenem__sucrose__6"),
                       c("biapenem__vehicle__3", "biapenem__glucose__3"),
                       c("biapenem__vehicle__6", "biapenem__glucose__6"),
                       c("biapenem__vehicle__3", "biapenem__fructose__3"),
                       c("biapenem__vehicle__6", "biapenem__fructose__6")),
    label = "p.signif", method = "wilcox.test", tip.length = 0.01, step.increase = 0.06) +
  labs(x = "Day", y = "Enterococcal\nCFU/gram") +
  theme_mouse() +
  legend_treatment()

save_panel(p_9c, "E9c_alternate_sugars.pdf", width = 9.5, height = 3.6)

# E9b: commercial smoothie (R45) ----------------------------------------------
groups_smoothie <- c("PBS_vehicle", "PBS_smoothie", "biapenem_vehicle", "biapenem_smoothie")
days_smoothie <- c(1, 3, 6)

smoothie <- read_mouse_sheet("Sup_Figure_9b_smoothie") |>
  rename(cfu = `CFUs per gram stool`, grp = Treatment, mouse = `Treatment + Mouse Number`) |>
  add_day("Treatement + Day") |>   # original sheet header is misspelled
  mutate(grp = factor(grp, levels = groups_smoothie),
         xvar = factor(str_glue("{grp}__{day}"), levels = xvar_levels(groups_smoothie, days_smoothie)))

p_9b <- smoothie |>
  ggplot(aes(x = xvar, y = cfu, colour = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_line(aes(group = mouse), colour = "gray80", linewidth = 0.4, alpha = 0.5) +
  geom_point(alpha = 0.7, shape = 16, size = 1) +
  scale_colour_manual(values = pal_smoothie4, name = "Treatment group", labels = pretty_grp) +
  scale_log10_sci() +
  scale_x_discrete(labels = rep(days_smoothie, length(groups_smoothie))) +
  stat_compare_means(
    comparisons = list(c("biapenem_vehicle__3", "biapenem_smoothie__3"),
                       c("biapenem_vehicle__6", "biapenem_smoothie__6")),
    label = "p.signif", method = "wilcox.test", tip.length = 0.04, step.increase = 0.15) +
  labs(x = "Day", y = "Enterococcal\nCFU/gram") +
  theme_mouse() +
  legend_treatment()

save_panel(p_9b, "E9b_smoothie.pdf", width = 5.2, height = 3.2)

# E9d: germ-free monocolonization time course ---------------------------------
mono_levels <- c("E. faecalis + Sucrose Water", "E. faecalis + Regular Water",
                 "PBS + Sucrose Water", "PBS + Regular Water")
pal_mono <- c("E. faecalis + Sucrose Water" = "#F8766D",
              "E. faecalis + Regular Water" = "#00BFC4",
              "PBS + Sucrose Water" = "#FAD5D2",
              "PBS + Regular Water" = "#A2E8E6")

mono <- read_mouse_sheet("Sup_Figure_9d_time_course_monoc") |>
  rename(hours = `Time post innoculation (hours)`, cfu = `CFUs per gram stool`,
         mouse = `Mouse Identifier`, treatment = Treatment) |>
  filter(!is.na(cfu)) |>
  mutate(treatment = factor(treatment, levels = mono_levels))

p_9d_tc <- mono |>
  ggplot(aes(x = hours, y = cfu, group = mouse, colour = treatment)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.2, alpha = 0.8, shape = 16) +
  scale_colour_manual(values = pal_mono) +
  scale_log10_sci() +
  scale_x_continuous(breaks = sort(unique(mono$hours))) +
  labs(x = "Hours post inoculation", y = "CFU/gram") +
  theme_mouse() +
  theme(legend.position = "right", legend.title = element_blank())

save_panel(p_9d_tc, "E9d_monocolonization_timecourse.pdf", width = 5.2, height = 3.2)

# E9d: AUC boxplot ------------------------------------------------------------
# Published comparison was a two-way ANOVA (ns); the Wilcoxon bracket here is the
# in-panel equivalent and will likewise read ns.
mono_auc <- read_mouse_sheet("Sup_Figure_9d_auc_boxplot") |>
  rename(auc = `Enterococcus load (AUC)`, treatment = Treatment) |>
  mutate(treatment = factor(treatment, levels = c("E. faecalis + Regular Water",
                                                  "E. faecalis + Sucrose Water")))

p_9d_auc <- mono_auc |>
  ggplot(aes(x = treatment, y = auc, colour = treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, shape = 16, size = 1) +
  scale_colour_manual(values = c("#00BFC4", "#F8766D")) +
  stat_compare_means(comparisons = list(c("E. faecalis + Regular Water",
                                          "E. faecalis + Sucrose Water")),
                     label = "p.signif", method = "wilcox.test", tip.length = 0.04) +
  labs(x = NULL, y = "Enterococcus load (AUC)") +
  theme_mouse() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_panel(p_9d_auc, "E9d_monocolonization_auc.pdf", width = 3.2, height = 3.4)

# E9h: bone-marrow-transplant recipients (R02) --------------------------------
# The published panel shows the BM-only (T-cell-depleted) recipients; the sheet
# also carries BM+Tcells, kept out of this panel to match the figure.
bmt_levels <- c("BMonly__PBS__vehicle", "BMonly__PBS__sucrose",
                "BMonly__biapenem__vehicle", "BMonly__biapenem__sucrose")
days_bmt <- c(0, 4, 6, 9)

bmt <- read_mouse_sheet("Sup_Figure_9h_bmt_cfu") |>
  rename(cfu = `CFUs per gram stool`, grp = `Treatment Group`) |>
  add_day("Treatment + Day") |>
  filter(str_starts(grp, "BMonly")) |>
  mutate(grp = factor(grp, levels = bmt_levels),
         xvar = factor(str_glue("{grp}__{day}"), levels = xvar_levels(bmt_levels, days_bmt)))

p_9h <- bmt |>
  ggplot(aes(x = xvar, y = cfu, colour = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, shape = 16, size = 0.9) +
  scale_colour_manual(values = pal_sucrose4, name = "Treatment group", labels = pretty_grp) +
  scale_log10_sci() +
  scale_x_discrete(labels = rep(days_bmt, length(bmt_levels))) +
  stat_compare_means(
    comparisons = list(c("BMonly__PBS__vehicle__4", "BMonly__PBS__sucrose__4"),
                       c("BMonly__PBS__vehicle__4", "BMonly__biapenem__vehicle__4"),
                       c("BMonly__biapenem__vehicle__6", "BMonly__biapenem__sucrose__6"),
                       c("BMonly__biapenem__vehicle__9", "BMonly__biapenem__sucrose__9")),
    label = "p.signif", method = "wilcox.test", tip.length = 0.02, step.increase = 0.12) +
  labs(x = "Day", y = "Enterococcal\nCFU/gram") +
  theme_mouse() +
  legend_treatment()

save_panel(p_9h, "E9h_bmt_cfu.pdf", width = 5.4, height = 3.4)
