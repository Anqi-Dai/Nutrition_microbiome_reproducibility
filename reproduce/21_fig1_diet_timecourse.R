# Figure 1 b,c,i,j,k,l,m: the cohort overview and diet/microbiome timecourses.
#
# Refactor of 072_Figure_1__code_for_Figure_1.Rmd (the histogram + loess + grid
# halves; the TaxUMAP panels e-h live in 20_fig1_taxumap.R). The l/m grids follow
# the corrected R64_F1_scatter.Rmd, which fixes a zero-fill bug in 072 (see the
# F1l/F1m block) and colours the food-group lines/strips. All panels are
# deterministic given the released tables, so nothing is cached. Panels:
#   F1b  meals recorded over transplant day (histogram)
#   F1c  stool samples over transplant day (histogram)
#   F1i  daily caloric intake over transplant day (loess)
#   F1j  diet alpha-diversity (Faith PD) over transplant day (loess)
#   F1k  microbiome alpha-diversity (inverse Simpson) over transplant day (loess)
#   F1l  per-day food-group dehydrated weight, 3x3 grid (loess)
#   F1m  per-day macronutrient grams, five-panel grid (loess)
#
# The published trend claims (each timecourse decreasing / increasing) are GEE
# Wald one-sided p-values; reproduced here and printed for verification.
#
# Inputs (all in released_data/): 152_combined_DTB.csv, 153_combined_META.csv,
# diet-alpha-diversity.tsv, food_group_color_key_final.csv,
# 072_total_patients_zero_eating_days_pid.csv.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggrastr)
  library(ggh4x)
  library(geepack)
  library(here)
})

# Format an integer with thousands separators for the panel annotations.
fmt <- function(n) format(n, big.mark = ",", trim = TRUE)

released <- function(f) file.path(Sys.getenv("NUTRITION_DATA", unset = here::here("released_data")), f)
results_dir <- here::here("results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Style params, carried over from the original.
axis_title_size <- 10
stip_txt_size <- 8
point_size <- 0.2
scatter_transparency <- 0.1
diet_line_color <- "#E41A1C"
stool_line_color <- "blue"
day0_line_color <- "gray40"
strip_color <- "gray91"
transplant_day <- 0
dayzeroline_size <- 1

day0_line <- function() geom_vline(xintercept = transplant_day, linetype = "dashed",
                                   linewidth = dayzeroline_size, color = day0_line_color)

# Inputs ----------------------------------------------------------------------
dtb  <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE)
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)
key  <- read_csv(released("food_group_color_key_final.csv"), col_types = "ccccc")
zerodaysall <- read_csv(released("072_total_patients_zero_eating_days_pid.csv"), show_col_types = FALSE)

# diet-alpha-diversity.tsv ships keyed by sampleid ("P{pid}d{fdrt}"); split it back
# into pid + fdrt so it plots and models like the original per-day faith table.
faith <- read_tsv(released("diet-alpha-diversity.tsv"), show_col_types = FALSE) %>%
  rename(sampleid = 1) %>%
  extract(sampleid, into = c("pid", "fdrt"), regex = "^(P\\d+)d(-?\\d+)$",
          convert = TRUE, remove = FALSE)

# F1b: meals recorded ---------------------------------------------------------
# Cohort counts derived from the tracker rather than hardcoded: one patient per
# pid, one meal per (pid, Meal, fdrt), one food item per row.
meals_tbl <- dtb %>% distinct(pid, Meal, fdrt)
diet_label <- str_glue("{fmt(n_distinct(dtb$pid))} patients,\n",
                       "{fmt(nrow(meals_tbl))} total meals\n",
                       "{fmt(nrow(dtb))} food items")

diet_hist <- meals_tbl %>%
  gghistogram(x = "fdrt", xlab = "Transplant day", ylab = "Meals recorded",
              color = diet_line_color, fill = diet_line_color, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 50, 20)) +
  day0_line() +
  annotate("text", x = 30, y = 600, hjust = 0, size = 2.5, label = diet_label) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 1 / 1.5, axis.title = element_text(size = axis_title_size))
ggsave(file.path(results_dir, "F1b_meals_hist.pdf"), diet_hist, width = 4, height = 3)

# F1c: stool samples ----------------------------------------------------------
# One patient per pid, one fecal sample per metadata row.
stool_label <- str_glue("{fmt(n_distinct(meta$pid))} patients,\n",
                        "{fmt(nrow(meta))} fecal samples")

stool_hist <- meta %>%
  ggplot(aes(x = sdrt)) +
  geom_histogram(alpha = 0.5, fill = stool_line_color, color = stool_line_color) +
  day0_line() +
  labs(x = "Transplant day", y = "Stool samples") +
  scale_x_continuous(breaks = seq(0, 50, 20)) +
  annotate("text", x = 30, y = 60, hjust = 0, size = 2.5, label = stool_label) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 1 / 1.5, axis.title = element_text(size = axis_title_size))
ggsave(file.path(results_dir, "F1c_stool_hist.pdf"), stool_hist, width = 4, height = 3)

# Shared loess timecourse: rasterised scatter + loess + day-0 line, sqrt y. -----
timecourse <- function(df, x, y, line_col, fill_col, title, xlab, ylab) {
  ggplot(df) +
    rasterise(geom_point(aes({{ x }}, {{ y }}), alpha = 0.3, size = point_size, shape = 16), dpi = 300) +
    geom_smooth(aes({{ x }}, {{ y }}), method = "loess", formula = "y ~ x",
                colour = line_col, linewidth = 1, fill = fill_col) +
    day0_line() +
    labs(x = xlab, y = ylab, title = title) +
    scale_x_continuous(breaks = seq(0, 50, 20)) +
    scale_y_sqrt() +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(size = axis_title_size),
          axis.title = element_text(size = axis_title_size), aspect.ratio = 1 / 1.5)
}

# F1i: caloric intake ---------------------------------------------------------
day_calori <- dtb %>%
  group_by(pid, fdrt) %>%
  summarise(daycal = sum(Calories_kcal), .groups = "drop")

cal_line <- day_calori %>%
  mutate(daycal = daycal / 1000) %>%
  timecourse(fdrt, daycal, diet_line_color, "hotpink", "Caloric intake", "", "*1000 Kcal")
ggsave(file.path(results_dir, "F1i_caloric_intake.pdf"), cal_line, width = 4, height = 3)

# F1j: diet alpha-diversity (Faith PD) ----------------------------------------
faith_line <- faith %>%
  mutate(faith_pd = faith_pd / 1000) %>%
  timecourse(fdrt, faith_pd, diet_line_color, "hotpink",
             "Diet a diversity", "Transplant day", "*1000")
ggsave(file.path(results_dir, "F1j_diet_diversity.pdf"), faith_line, width = 4, height = 3)

# F1k: microbiome alpha-diversity (inverse Simpson) ---------------------------
stool_alpha <- meta %>%
  timecourse(sdrt, simpson_reciprocal, stool_line_color, "darkblue",
             "Microbiome\nalpha diversity", "", "")
ggsave(file.path(results_dir, "F1k_microbiome_diversity.pdf"), stool_alpha, width = 4, height = 3)

# Per-day grams, kept long. Zeros are added ONLY for the confirmed zero-eating
# days (cross-joined with the group names), not for every group missing from an
# eating day: the earlier spread(fill = 0) injected a 0 for each food group not
# eaten on each day, which dragged the loess down (the R64 correction).
zero_days <- zerodaysall %>% select(pid, fdrt)

# F1m: macronutrient grid -----------------------------------------------------
m_all <- dtb %>%
  select(pid, fdrt, Protein_g:Sugars_g) %>%
  pivot_longer(Protein_g:Sugars_g, names_to = "grp", values_to = "gram") %>%
  mutate(grp = str_replace(grp, "_g$", "")) %>%
  group_by(pid, fdrt, grp) %>%
  summarise(eachsum = sum(gram), .groups = "drop")

dailymacro <- bind_rows(
  m_all,
  zero_days %>% cross_join(tibble(grp = unique(m_all$grp))) %>% mutate(eachsum = 0))

m_panel <- dailymacro %>%
  mutate(grp = if_else(str_detect(grp, "Carbohydrates"), "Carbs", grp)) %>%
  ggplot() +
  rasterise(geom_point(aes(fdrt, eachsum), alpha = scatter_transparency, size = point_size, shape = 16), dpi = 300) +
  geom_smooth(aes(fdrt, eachsum), method = "loess", formula = "y ~ x",
              colour = diet_line_color, linewidth = 1, fill = "hotpink") +
  day0_line() +
  labs(x = "Transplant day", y = "Grams", title = "Macronutrients") +
  facet_wrap(~ grp, nrow = 3, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 50, 20)) +
  scale_y_sqrt() +
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        strip.background = element_rect(color = "white", fill = strip_color, linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = stip_txt_size),
        axis.title = element_text(size = axis_title_size),
        plot.title = element_text(size = axis_title_size), aspect.ratio = 1)
ggsave(file.path(results_dir, "F1m_macronutrients.pdf"), m_panel, width = 5, height = 6)

# F1l: food-group grid --------------------------------------------------------
fg_all <- dtb %>%
  mutate(fgrp1 = str_sub(as.character(Food_code), 1, 1)) %>%
  group_by(pid, fdrt, fgrp1) %>%
  summarise(eachsum = sum(dehydrated_weight), .groups = "drop") %>%
  left_join(key %>% select(fgrp1, shortname), by = "fgrp1") %>%
  select(pid, fdrt, grp = shortname, eachsum)

fg_total <- bind_rows(
  fg_all,
  zero_days %>% cross_join(tibble(grp = unique(fg_all$grp))) %>% mutate(eachsum = 0))

# Each food group's loess line + ribbon takes its own colour; the facet strips are
# painted the same colour (alphabetical, to match the facet layout) via ggh4x.
fg_colors <- setNames(key$color, key$shortname)
strip_colors <- key %>%
  filter(shortname %in% unique(fg_total$grp)) %>%
  arrange(shortname) %>%
  pull(color)

fg_panel <- fg_total %>%
  ggplot() +
  rasterise(geom_point(aes(fdrt, eachsum), alpha = scatter_transparency, size = point_size, shape = 16), dpi = 300) +
  geom_smooth(aes(fdrt, eachsum, color = grp, fill = grp), method = "loess",
              formula = "y ~ x", linewidth = 1) +
  day0_line() +
  labs(x = "Transplant day", y = "Grams", title = "Food groups") +
  facet_wrap2(~ grp, scales = "free_y", nrow = 3,
              strip = strip_themed(background_x = elem_list_rect(fill = strip_colors))) +
  scale_color_manual(values = fg_colors) +
  scale_fill_manual(values = fg_colors) +
  scale_x_continuous(breaks = seq(0, 50, 20)) +
  scale_y_sqrt() +
  theme_classic() +
  theme(axis.text = element_text(size = 6),
        strip.text.x = element_text(size = stip_txt_size),
        axis.title = element_text(size = axis_title_size),
        plot.title = element_text(size = axis_title_size),
        aspect.ratio = 1, legend.position = "none")
ggsave(file.path(results_dir, "F1l_food_groups.pdf"), fg_panel, width = 6, height = 6)

# Trend p-values (GEE, AR1 working corr, one-sided Wald) -----------------------
# id = pid must resolve as a symbol inside the data frame (passing it as an
# external vector makes geepack's AR1 estimation pathologically slow).
gee_pval <- function(df, timevar) {
  d <- df %>% mutate(pid = as.factor(pid))
  mod <- geeglm(reformulate(timevar, "y"), family = gaussian, corstr = "ar1",
                id = pid, data = d)
  z <- coef(mod)[2] / sqrt(vcov(mod)[2, 2])
  1 - pnorm(abs(z))
}

p_cal   <- gee_pval(day_calori %>% mutate(y = log(daycal + 1)), "fdrt")
p_faith <- gee_pval(faith %>% mutate(y = log(faith_pd)), "fdrt")
p_stool <- gee_pval(meta %>% mutate(y = log(simpson_reciprocal)), "sdrt")

# Per-group macro / food-group trends, restricted to the early window (fdrt <= 12).
grp_pvals <- function(long_df) {
  long_df %>%
    filter(fdrt <= 12) %>%
    group_split(grp) %>%
    map_dfr(function(d) tibble(grp = d$grp[1],
                               pval = gee_pval(d %>% mutate(y = log(eachsum + 0.5)), "fdrt"))) %>%
    mutate(FDR = p.adjust(pval, method = "BH"), q_lt_005 = FDR < 0.05)
}
macro5 <- grp_pvals(dailymacro)
fg9    <- grp_pvals(fg_total)

message("F1 b,c,i,j,k,l,m written to results/.")
message(sprintf("GEE one-sided p: caloric=%.3g  diet-faith=%.3g  microbiome=%.3g",
                p_cal, p_faith, p_stool))
message("Macronutrient trend p-values (fdrt<=12):")
print(macro5)
message("Food-group trend p-values (fdrt<=12):")
print(fg9)
