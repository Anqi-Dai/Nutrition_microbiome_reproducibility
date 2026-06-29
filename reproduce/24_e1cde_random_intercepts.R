# Extended Fig. E1 c,d,e: patient-specific diversity deviations (random intercepts
# of the F2d food-group diversity model) and their relationship to sweets intake.
#
# Refactor of R53. The random intercepts r_pid[<pid>,Intercept] are the per-patient
# "baseline" diversity deviations after the fixed effects (antibiotics, diet,
# conditioning intensity) are subtracted. We rank them, and ask whether average
# sweets intake still tracks them (leftover sweet-related variance leaking into the
# random effects).
#
#   E1c  ranked forest of patient random intercepts (median + 95% CrI), uncoloured
#   E1d  scatter of average daily sweets intake vs the median random intercept,
#        with a Spearman correlation
#   E1e  the same forest, coloured by log2(avg daily sweets + 1)
#
# Input: the posterior draws of the F2d main model. R53 read a precomputed
# 172_div_fg_posterior_draws.csv; here that file is generated once from the cached
# F2d fit (intermediate_data/172_fit_fg_diversity.rds, from 10_fit_diversity_models.R)
# -- it IS the F2d main model -- and cached to intermediate_data/ for reuse.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages(library(posterior))

dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE)

# ---- 1. posterior random-intercept draws from the F2d fit ------------------
draws_csv <- cache_path("172_div_fg_posterior_draws.csv")
if (!file.exists(draws_csv)) {
  message("extracting r_pid draws from the cached F2d fit ...")
  fit <- readRDS(cache_path("172_fit_fg_diversity.rds"))
  rp <- as_draws_df(fit)
  rp <- rp[, grep("^r_pid\\[", variables(rp))]
  readr::write_csv(as_tibble(rp), draws_csv)
}
posterior_draws_df <- read_csv(draws_csv, show_col_types = FALSE)

intercept_medians <- posterior_draws_df |>
  pivot_longer(starts_with("r_pid"), names_to = "parameter", values_to = "draw") |>
  mutate(pid = str_extract(parameter, "(?<=r_pid\\[)[^,]+")) |>
  filter(!is.na(pid)) |>
  group_by(pid) |>
  summarise(median_intercept = median(draw),
            lower_ci = quantile(draw, 0.025),
            upper_ci = quantile(draw, 0.975), .groups = "drop")

# ---- 2. average daily sweets intake per patient (food group 9) -------------
patient_avg_sweets <- dtb |>
  mutate(Food_code = as.character(Food_code)) |>
  group_by(pid) |>
  summarise(total_sweets_grams = sum(dehydrated_weight[str_starts(Food_code, "9")], na.rm = TRUE),
            total_eating_days = n_distinct(fdrt), .groups = "drop") |>
  mutate(avg_daily_sweets = total_sweets_grams / total_eating_days)

combined_data <- inner_join(patient_avg_sweets, intercept_medians, by = "pid") |>
  mutate(log2_avg_daily_sweets = log2(avg_daily_sweets + 1))

# ---- 3. Spearman correlation (E1d annotation) ------------------------------
ct <- cor.test(combined_data$avg_daily_sweets, combined_data$median_intercept,
               method = "spearman")
rho <- unname(ct$estimate); pval <- ct$p.value
p_string <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
annotation <- sprintf("Spearman's rho = %.2f\n%s", rho, p_string)
message(sprintf("E1d Spearman: rho = %.3f, %s (n = %d patients)",
                rho, p_string, nrow(combined_data)))

# ---- 4. panels (styled like the published E1 c,d,e) ------------------------
# Panel letters: c = uncoloured forest, d = forest coloured by sweets, e = scatter.
forest_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 12))
x_dev <- scale_x_continuous(breaks = seq(-1, 1, 0.5))
xlab_dev <- "Patient-specific deviation from\naverage microbiome ln(diversity)"
patient_order <- aes(x = median_intercept, y = reorder(pid, median_intercept),
                     xmin = lower_ci, xmax = upper_ci)

# c: uncoloured forest
forest_c <- ggplot(combined_data, patient_order) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(height = 0, linewidth = 0.4) +
  geom_point(size = 0.7) +
  x_dev + labs(x = xlab_dev, y = "individual patients") + forest_theme

# d: forest coloured by log2 average daily sweets
forest_d <- ggplot(combined_data,
    aes(x = median_intercept, y = reorder(pid, median_intercept),
        xmin = lower_ci, xmax = upper_ci, color = log2_avg_daily_sweets)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(height = 0, linewidth = 0.5) +
  geom_point(size = 0.7, color = "grey25") +
  scale_color_viridis_c(option = "turbo", breaks = c(2, 4, 6),
                        labels = c("2  low", "4", "6  high"),
                        name = expression(atop("sweets", "intake (log"[2] * ")"))) +
  x_dev + labs(x = xlab_dev, y = "individual patients") + forest_theme +
  theme(legend.position = c(0.85, 0.25), legend.title = element_text(size = 10))

# e: scatter of sweets intake vs patient intercept (Spearman)
ann <- bquote(atop("Spearman's" ~ rho == .(sprintf("%.2f", rho)),
                   "p" == .(sprintf("%.2f", pval))))
scatter_e <- ggplot(combined_data, aes(avg_daily_sweets, median_intercept)) +
  geom_point(alpha = 0.7, size = 2, color = "grey20") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", formula = y ~ x) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, size = 4.2,
           color = "#9A7D0A", parse = TRUE,
           label = as.character(as.expression(ann))) +
  labs(x = "Daily Sweets Intake (averaged g/day)",
       y = "Patient-Specific Intercept\n(posterior median)") +
  theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())

# ---- 5. save each panel independently --------------------------------------
out_dir <- here::here("results"); if (!dir.exists(out_dir)) dir.create(out_dir)
ggsave(file.path(out_dir, "E1c_random_intercepts.pdf"),        forest_c, width = 5, height = 7, units = "in")
ggsave(file.path(out_dir, "E1d_random_intercepts_sweets.pdf"), forest_d, width = 5.4, height = 7, units = "in")
ggsave(file.path(out_dir, "E1e_sweets_correlation.pdf"),       scatter_e, width = 6.5, height = 4.5, units = "in")
message("wrote results/E1c_random_intercepts.pdf, E1d_random_intercepts_sweets.pdf, ",
        "E1e_sweets_correlation.pdf")
