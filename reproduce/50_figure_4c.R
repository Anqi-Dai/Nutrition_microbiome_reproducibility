# Figure 4c. Sucrose exacerbates Enterococcus expansion after antibiotic.
#
# Two stacked sub-panels share the same four treatment groups (PBS vs biapenem,
# crossed with vehicle vs sucrose):
#   F4c upper  fecal enterococcal burden (CFU/g) on days 1, 3, 6
#   F4c lower  the same burden summarised per mouse as a trapezoidal AUC
#
# Style follows the originals (R07 / R45): grey-vs-pink palette, log10 axis with
# 10^n ticks, days relabelled under each treatment block, and the key Wilcoxon
# comparisons drawn in-panel with stat_compare_means. This is the template the
# other CFU scripts follow.

source(here::here("reproduce", "mouse", "_mouse_helpers.R"))

groups_4c <- c("PBS__vehicle", "PBS__sucrose", "biapenem__vehicle", "biapenem__sucrose")

# F4c upper: CFU over days ----------------------------------------------------
indiv <- read_mouse_sheet("Figure_4c_indiv_days_data") |>
  rename(cfu = `CFUs per gram stool`, grp = `Treatment Group`) |>
  add_day("Treatment + Day") |>
  mutate(grp = factor(grp, levels = groups_4c),
         xvar = factor(str_glue("{grp}__{day}"),
                       levels = xvar_levels(groups_4c, sort(unique(day)))))

p_indiv <- indiv |>
  ggplot(aes(x = xvar, y = cfu, colour = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, shape = 16, size = 1) +
  scale_colour_manual(values = pal_sucrose4, name = "Treatment group", labels = pretty_grp) +
  scale_log10_sci() +
  scale_x_discrete(labels = rep(sort(unique(indiv$day)), length(groups_4c))) +
  stat_compare_means(
    comparisons = list(c("biapenem__vehicle__3", "biapenem__sucrose__3"),
                       c("biapenem__vehicle__6", "biapenem__sucrose__6")),
    label = "p.signif", method = "wilcox.test", tip.length = 0.02,
    step.increase = 0.12) +
  labs(x = "Day", y = "Enterococcal\nCFU/gram") +
  theme_mouse() +
  legend_treatment()

save_panel(p_indiv, "F4c_cfu_over_days.pdf", width = 6.0, height = 3.2)

# F4c lower: trapezoidal AUC --------------------------------------------------
auc <- read_mouse_sheet("Figure_4c_trapezoidal_auc") |>
  rename(auc = `Trapezoidal AUC Value`, grp = `Treatment Group`) |>
  mutate(grp = factor(grp, levels = groups_4c))

p_auc <- auc |>
  ggplot(aes(x = grp, y = auc, colour = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5, shape = 16, size = 1) +
  scale_colour_manual(values = pal_sucrose4) +
  scale_log10_sci() +
  stat_compare_means(
    comparisons = list(c("PBS__vehicle", "biapenem__vehicle"),
                       c("biapenem__vehicle", "biapenem__sucrose"),
                       c("PBS__vehicle", "PBS__sucrose")),
    label = "p.signif", method = "wilcox.test", tip.length = 0.04) +
  labs(x = NULL, y = "Trapezoidal\nAUC") +
  theme_mouse() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_panel(p_auc, "F4c_auc.pdf", width = 3.2, height = 3.4)
