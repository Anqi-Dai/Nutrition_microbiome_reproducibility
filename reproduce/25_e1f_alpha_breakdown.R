# Extended Fig. E1f: unadjusted alpha-diversity (ln inverse Simpson) broken down by
# conditioning intensity, by antibiotics x sweets, and by TPN x sweets.
#
# Refactor of R44 (data logic kept exactly), restyled toward the published panel:
# grey nested section strips on the left (the antibiotics block split into
# not-exposed / exposed sub-strips, as in the published bracket), black boxplots
# over grey jittered points, an "n samples" column on the right, and
# x = ln(inverse Simpson index).
#
# Groupings (all from 153_combined_META):
#   - conditioning intensity: ablative -> "myeloablative", reduced, nonablative
#   - antibiotics (empirical) x sweets: sweets split at the median intake AMONG
#     consumers (fg_sweets > 0); no / below-median / above-median, within
#     not-antibiotic-exposed and antibiotic-exposed sub-strips
#   - TPN x any sweets, excluding EN-exposed samples (EN counts as sweets)
# Verified the group sizes match the published panel exactly (629/288/92;
# 41/276/253/29/194/216; 35/8/893/62).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggh4x)   # facet_nested: nested strips + free/proportional panel heights
})

meta <- read_csv(here::here("released_data", "153_combined_META.csv"), show_col_types = FALSE) |>
  mutate(value = log(simpson_reciprocal))
med_sweets <- meta |> filter(fg_sweets > 0) |> summarise(m = median(fg_sweets)) |> pull(m)

sec <- c(int = "Unadjusted α-diversity by\nconditioning intensity",
         abx = "by exposure to\nantibiotics and sweets",
         tpn = "by exposure to total parenteral\nnutrition (TPN) and sweets")
sweet_cat <- function(fg) case_when(fg == 0 ~ "no sweets",
                                    fg <= med_sweets ~ "below-median sweets",
                                    TRUE ~ "above-median sweets")

d_int <- meta |> transmute(value, section = sec["int"], sub = " ",
  row = recode(intensity, ablative = "myeloablative", reduced = "reduced",
               nonablative = "nonablative"))
d_abx <- meta |> transmute(value, section = sec["abx"], row = sweet_cat(fg_sweets),
  sub = if_else(empirical, "antibiotic-\nexposed", "not antibiotic-\nexposed"))
d_tpn <- meta |> filter(!EN) |> transmute(value, section = sec["tpn"], sub = "  ",
  row = case_when(TPN & fg_sweets > 0 ~ "TPN, Sweets", TPN ~ "TPN, No Sweets",
                  fg_sweets > 0 ~ "No TPN, Sweets", TRUE ~ "No TPN, No Sweets"))

row_levels <- c("nonablative", "reduced", "myeloablative",
                "above-median sweets", "below-median sweets", "no sweets",
                "No TPN, No Sweets", "No TPN, Sweets", "TPN, No Sweets", "TPN, Sweets")
plot_df <- bind_rows(d_int, d_abx, d_tpn) |>
  mutate(section = factor(section, levels = sec),
         sub = factor(sub, levels = c(" ", "not antibiotic-\nexposed",
                                      "antibiotic-\nexposed", "  ")),
         row = factor(row, levels = row_levels))

n_df <- plot_df |> count(section, sub, row, name = "n_samples")
n_x <- max(plot_df$value) + 0.55                       # n-samples column position

f <- ggplot(plot_df, aes(value, row)) +
  geom_jitter(height = 0.28, width = 0, size = 0.5, alpha = 0.18, colour = "grey30") +
  geom_boxplot(fill = "white", colour = "black", linewidth = 0.45,
               width = 0.5, outlier.shape = NA, alpha = 0.6) +
  geom_text(data = n_df, aes(x = n_x, y = row, label = n_samples),
            colour = "grey45", size = 3.2) +
  # "n samples" header, nudged above the top (intensity) panel
  geom_text(data = tibble(section = factor(sec["int"], levels = sec), sub = factor(" "),
                          row = factor("myeloablative", levels = row_levels)),
            aes(x = n_x, y = row, label = "n samples"), vjust = -2.6,
            fontface = "bold", size = 3.3, inherit.aes = FALSE) +
  facet_nested(rows = vars(section, sub), scales = "free_y", space = "free_y",
               switch = "y") +
  scale_x_continuous(breaks = 0:3, limits = c(NA, n_x + 0.3),
                     expand = expansion(mult = c(0.01, 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "ln(inverse Simpson index)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 8.5),
        strip.placement = "outside",
        panel.spacing = unit(4, "mm"),
        plot.margin = margin(24, 10, 6, 6),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

out <- here::here("results", "E1f_alpha_diversity_breakdown.pdf")
ggsave(out, f, width = 9, height = 8, units = "in", device = cairo_pdf)
message("wrote ", out, "  (median sweets among consumers = ", round(med_sweets, 2), " g)")
