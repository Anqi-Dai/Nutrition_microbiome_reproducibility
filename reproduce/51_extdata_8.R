# Extended Data Figure E8, non-sequencing panels.
#
# E8 a-f come from 16S and are produced by 53_extdata_8_16s.R. This file covers
# the panels that read a self-contained source-data sheet and need no count table:
#   E8g  daily chow consumption per day (biapenem + sucrose vs + vehicle)
#   E8h  trapezoidal AUC of chow consumption
#   E8j  AUC of Enterococcus under delayed sugar, both experiments pooled
#   E8k  delayed-sucrose CFU/g time course, one facet per treatment
#   E8l  Enterococcus CFU/g on fibre-free chow, already in log10 units
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
  rename(auc = `AUC of CFU/gram *days`, grp = `Treatment Group`,
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
  labs(x = NULL, y = "Area Under the Curve\n(log10 CFU/gram * days)") +
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
  mutate(cfu = as.numeric(cf_us_per_gram_stool),  # column ships as text with literal "NA"
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
  labs(x = "Day of experiment", y = "Enterococcal\nCFU/gram") +
  theme_mouse() +
  theme(strip.text = element_text(size = 7))

save_panel(p_8k, "E8k_cfu_over_days.pdf", width = 8.0, height = 3.0)

# E8l: fibre-free chow, log10 CFU (R07) ---------------------------------------
groups_dpbs <- c("DPBS__vehicle", "DPBS__sucrose", "biapenem__vehicle", "biapenem__sucrose")
days_8l <- c(1, 3, 6, 9)

cfu_8l <- read_mouse_sheet("Sup_Figure_8l_no_fiber_chow_CFU") |>
  rename(log_cfu = `log 10 (Enterococcal CFU/gram stool)`, grp = `Treatment Group`) |>
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
  labs(x = "Day", y = expression(log[10] ~ "Enterococcal CFU/gram")) +
  theme_mouse() +
  legend_treatment()

save_panel(p_8l, "E8l_fiberfree_log10_cfu.pdf", width = 6.0, height = 3.2)

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
