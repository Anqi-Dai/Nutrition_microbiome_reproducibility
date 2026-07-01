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
# Group-major / day-minor layout (R05 f3_five_days): the four antibiotic x diet
# arms sit as blocks, days ascending within each, coloured by group. The brackets
# are the R05 Wilcoxon contrasts: biapenem vehicle vs biapenem sucrose, and
# biapenem vehicle vs PBS vehicle, at each of days 3, 6, 9, 14, 21.
grp_9a <- c("PBS__vehicle", "PBS__sucrose", "biapenem__vehicle", "biapenem__sucrose")
days_9a <- c(1, 3, 6, 9, 14, 21)

course_9a <- read_course("Sup_Figure_9a_21_day_exp") |>
  mutate(antibiotic = recode(antibiotic, DPBS = "PBS"),
         day = as.integer(as.character(day)),
         grp = factor(str_glue("{antibiotic}__{sugar}"), levels = grp_9a),
         xvar = factor(str_glue("{antibiotic}__{sugar}__{day}"),
                       levels = xvar_levels(grp_9a, days_9a)))

cmp_9a <- c(
  lapply(days_9a[-1], function(d) c(str_glue("biapenem__vehicle__{d}"), str_glue("biapenem__sucrose__{d}"))),
  lapply(days_9a[-1], function(d) c(str_glue("biapenem__vehicle__{d}"), str_glue("PBS__vehicle__{d}")))
)

p_9a <- course_9a |>
  ggplot(aes(x = xvar, y = cfu, colour = grp)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, shape = 16, size = 0.8) +
  scale_colour_manual(values = setNames(pal_sucrose4, grp_9a),
                      name = "Antibiotic + diet", labels = pretty_grp) +
  scale_log10_sci() +
  scale_x_discrete(labels = rep(days_9a, length(grp_9a))) +
  stat_compare_means(comparisons = cmp_9a, label = "p.signif", method = "wilcox.test",
                     method.args = list(exact = TRUE, correct = TRUE),
                     tip.length = 0.01, step.increase = 0.06, vjust = 0.6) +
  labs(x = "Day", y = "Enterococcal\nCFU/gram") +
  theme_mouse() +
  legend_treatment()

save_panel(p_9a, "E9a_21day_timecourse.pdf", width = 9.5, height = 3.6)

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
# Faceted by water (Sucrose / Regular), with E. faecalis (red) vs PBS (black) as
# the two colours. Each mouse is a faint thin line; the thick line is the
# per-hour median. CFU is on a log axis; the PBS arm sits at the 10^0 floor,
# labelled "0" to read as below-detection.
mono_water_levels <- c("Sucrose Water", "Regular Water")
pal_microbe <- c("E. faecalis" = "red", "PBS" = "black")

mono <- read_mouse_sheet("Sup_Figure_9d_time_course_monoc") |>
  rename(hours = `Time post innoculation (hours)`, cfu = `CFUs per gram stool`,
         mouse = `Mouse Identifier`, treatment = Treatment) |>
  filter(!is.na(cfu)) |>
  separate(treatment, into = c("microbe", "water"), sep = " \\+ ", remove = FALSE) |>
  mutate(water = factor(water, levels = mono_water_levels),
         microbe = factor(microbe, levels = c("E. faecalis", "PBS")))

mono_med <- mono |>
  group_by(water, microbe, hours) |>
  summarise(cfu = median(cfu), .groups = "drop")

p_9d_tc <- ggplot(mono, aes(x = hours, y = cfu, colour = microbe)) +
  geom_line(aes(group = mouse), alpha = 0.25, linewidth = 0.4) +
  geom_point(alpha = 0.3, size = 1, shape = 16) +
  geom_line(data = mono_med, aes(group = microbe), linewidth = 1.3) +
  facet_wrap(~water) +
  scale_colour_manual(values = pal_microbe, name = NULL) +
  scale_y_log10(breaks = c(1, 1e4, 1e8, 1e12),
                labels = c(expression(0), expression(10^4), expression(10^8), expression(10^12))) +
  scale_x_continuous(breaks = c(0, 4, 8, 24, 48, 144)) +
  labs(x = "hours since inoculation", y = "Enterococcus\nCFU / g feces") +
  theme_mouse() +
  theme(legend.position = c(0.78, 0.28),
        legend.text = element_text(size = 10),
        legend.key = element_blank(), legend.background = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 3)))

save_panel(p_9d_tc, "E9d_monocolonization_timecourse.pdf", width = 6.6, height = 3.2)

# E9d: AUC boxplot ------------------------------------------------------------
# AUC is computed per mouse from the CFU time course exactly as R33: area under
# the log10(CFU)-vs-hours curve (DescTools::AUC) over the E. faecalis arms.
# experiment_no is recovered from the E1/E2 mouse-id prefix so R33's two-way
# ANOVA (Treatment + experiment_no, p = 0.7) can be run. Load is shown as
# AUC / 1000 on a single y axis that starts at 0, per the panel.
mono_auc <- mono |>
  filter(microbe == "E. faecalis") |>
  mutate(experiment_no = factor(str_extract(mouse, "^E[0-9]")),
         water = factor(recode(as.character(water),
                               "Sucrose Water" = "sucrose\nwater",
                               "Regular Water" = "regular\nwater"),
                        levels = c("sucrose\nwater", "regular\nwater"))) |>
  group_by(mouse, water, experiment_no) |>
  summarise(auc = DescTools::AUC(x = hours, y = log10(cfu)), .groups = "drop") |>
  mutate(load = auc / 1000)

auc_p <- summary(aov(auc ~ water + experiment_no, data = mono_auc))[[1]][["Pr(>F)"]][1]

p_9d_auc <- mono_auc |>
  ggplot(aes(x = water, y = load)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "#f6a6a0", colour = "gray25") +
  geom_jitter(width = 0.15, alpha = 0.7, shape = 16, size = 1.6, colour = "gray25") +
  annotate("segment", x = 1, xend = 2, y = 1.55, yend = 1.55, linewidth = 0.5) +
  annotate("text", x = 1.5, y = 1.58, label = paste0("p = ", signif(auc_p, 1)),
           size = 4, vjust = 0) +
  scale_y_continuous(limits = c(0, 1.65), expand = c(0, 0)) +
  labs(x = NULL, y = expression(atop("Enterococcus load", "(AUC x" * 10^3 * ")"))) +
  theme_mouse()

save_panel(p_9d_auc, "E9d_monocolonization_auc.pdf", width = 3.0, height = 3.4)

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
