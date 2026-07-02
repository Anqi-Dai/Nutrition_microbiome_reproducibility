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
#   Supplementary Tables 1-6 in one workbook (S1 full cohort, S2 microbiome
#        sub-cohort, S3/S5 characteristics by sugar / cluster, S4/S6 Cox models)
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
  library(openxlsx)
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

# Assemble one workbook, one sheet per table: short title, table, caption.
tables <- list(
  list(tbl = s1, sheet = "Supplementary Table 1",
       title = "Supplementary Table 1. Patient characteristics of the full nutrition cohort.",
       caption = "Data are presented as mean (standard deviation) for continuous variables and as count (percentage) for categorical variables. PBSC, peripheral blood stem cells. MDS/MPN, myelodysplastic syndromes/myeloproliferative neoplasms. GVHD, graft-versus-host disease. CNI, calcineurin inhibitor. PTCy, post-transplant cyclophosphamide."),
  list(tbl = s2, sheet = "Supplementary Table 2",
       title = "Supplementary Table 2. Patient characteristics of the microbiome-nutrition sub-cohort.",
       caption = "Data are presented as mean (standard deviation) for continuous variables and as count (percentage) for categorical variables. PBSC, peripheral blood stem cells. MDS/MPN, myelodysplastic syndromes/myeloproliferative neoplasms. GVHD, graft-versus-host disease. CNI, calcineurin inhibitor. PTCy, post-transplant cyclophosphamide."),
  list(tbl = s3, sheet = "Supplementary Table 3",
       title = "Supplementary Table 3. Patient characteristics by dichotomized averaged sugar intake, normalized by total caloric intake.",
       caption = "Patients dichotomized at the median sugar intake (above- vs below-median), normalized by total caloric intake and averaged across HCT days -7 to 12. Median (IQR) for continuous, count (percentage) for categorical variables."),
  list(tbl = s4, sheet = "Supplementary Table 4",
       title = "Supplementary Table 4. Multivariable Cox proportional hazards model for overall survival (OS) in dichotomized sugar-intake groups.",
       caption = "Adjusted HR and 95% CI, landmarked at day 12 after HCT, with an interaction between sugar-intake group and days of broad-spectrum antibiotic exposure (day -7 to 12), adjusting for conditioning intensity, graft type, and GVHD prophylaxis. The Abx-exposure row is the HR in the above-median group."),
  list(tbl = s5, sheet = "Supplementary Table 5",
       title = "Supplementary Table 5. Patient characteristics by diet-pattern cluster.",
       caption = "Clusters were defined by latent-trajectory class analysis without clinical variables. Median (IQR) for continuous, count (percentage) for categorical variables."),
  list(tbl = s6, sheet = "Supplementary Table 6",
       title = "Supplementary Table 6. Multivariable Cox proportional hazards model for overall survival (OS) in dietary-pattern clusters.",
       caption = "Adjusted HR and 95% CI, landmarked at day 12 after HCT, with an interaction between diet-pattern cluster and days of broad-spectrum antibiotic exposure, adjusting for conditioning intensity, graft type, and GVHD prophylaxis. The Abx-exposure row is the HR in Cluster 1.")
)

wb <- createWorkbook()
title_style  <- createStyle(textDecoration = "bold", fontSize = 12)
header_style <- createStyle(textDecoration = "bold", fgFill = "#D9D9D9", border = "bottom")
cap_style    <- createStyle(fontSize = 9, textDecoration = "italic")

for (t in tables) {
  addWorksheet(wb, t$sheet)
  writeData(wb, t$sheet, t$title, startRow = 1)
  addStyle(wb, t$sheet, title_style, rows = 1, cols = 1)

  df <- t$tbl |> as_tibble()
  names(df) <- names(df) |> str_remove_all("\\*\\*|__")
  df <- df |> mutate(across(everything(), ~ str_remove_all(as.character(.x), "\\*\\*|__")))
  writeData(wb, t$sheet, df, startRow = 3)
  addStyle(wb, t$sheet, header_style, rows = 3, cols = seq_len(ncol(df)), gridExpand = TRUE)

  cap_row <- nrow(df) + 5
  writeData(wb, t$sheet, t$caption, startRow = cap_row)
  addStyle(wb, t$sheet, cap_style, rows = cap_row, cols = 1)
  setColWidths(wb, t$sheet, cols = seq_len(ncol(df)), widths = "auto")
}
saveWorkbook(wb, file.path(results_dir, "Supplementary_Tables_1_6.xlsx"), overwrite = TRUE)

message("Wrote Fig 3 / E6 panels and Supplementary_Tables_1_6.xlsx to results/.")
