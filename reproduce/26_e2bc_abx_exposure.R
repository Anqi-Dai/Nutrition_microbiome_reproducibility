# Extended Fig. E2 b,c: sample-level antibiotic exposures (the two-day window prior
# to each fecal sample) over transplant day, as an absolute-count stacked histogram
# (E2b) and a 100% stacked percentage (E2c).
#
# Refactor of R21. Each fecal sample is assigned a single exposure category from its
# prior-2-day medication exposures (Data_S4), by the hierarchy broad-spectrum >
# fluoroquinolones > other antibacterials > not antibacterial. Counts/percentages
# are then tallied per transplant day (sdrt).
#
# E2a (the per-patient antibiotic-usage heatmap) is NOT reproducible from the
# released data -- it needs each patient's full daily antibiotic time course, which
# is not released; Data_S4 only carries the prior-2-day exposure per fecal sample.

suppressPackageStartupMessages(library(tidyverse))

released <- function(f) here::here("released_data", f)
out_dir  <- here::here("results"); if (!dir.exists(out_dir)) dir.create(out_dir)

expodat <- read_csv(released("Data_S4_Medication_Exposures_in_the_Two_Days_Prior_to_Stool_Sample_Collection.csv"),
                    show_col_types = FALSE)

# one exposure category per sample (hierarchy)
meds_categorized <- expodat |> distinct() |>
  group_by(sampleid, sdrt) |>
  summarise(exposure_category = case_when(
      any(drug_category_for_this_study == "broad_spectrum") ~ "broad_spectrum",
      any(drug_category_for_this_study == "fluoroquinolones") ~ "fluoroquinolones",
      any(drug_category_for_this_study == "other_antibacterials") ~ "other_antibacterials",
      TRUE ~ "not_antibacterial"), .groups = "drop")

meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)
total_by_day <- meta |> group_by(sdrt) |> summarise(daytotal = n(), .groups = "drop")

exposure_grps <- meds_categorized |>
  distinct(sampleid, exposure_category, sdrt) |>
  count(sdrt, exposure_category) |>
  inner_join(total_by_day, by = "sdrt") |>
  mutate(perc = round(n / daytotal * 100, 2),
         exposure_category = factor(exposure_category,
           levels = c("broad_spectrum", "fluoroquinolones",
                      "other_antibacterials", "not_antibacterial")))

# palette + legend labels matched to the published panel (soft salmon / periwinkle /
# pale pink / grey; R21 used tomato/steelblue/goldenrod3/gray as working colours).
pal <- c(broad_spectrum = "#E07B73", fluoroquinolones = "#A7AFD4",
         other_antibacterials = "#F2C6D8", not_antibacterial = "#BEBEBE")
lab <- c(broad_spectrum = "broad-spectrum",
         fluoroquinolones = "prophylactic fluoroquinolone",
         other_antibacterials = "other antibacterial",
         not_antibacterial = "no antibacterial exposure")

base <- list(
  scale_fill_manual(values = pal, labels = lab, name = "Antibiotic Exposure Category"),
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1),
  guides(fill = guide_legend(ncol = 1, title.position = "top")),
  labs(x = "Transplant day"),
  theme_minimal(base_size = 12),
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        legend.title = element_text(face = "bold")))

# E2b: absolute counts
e2b <- ggplot(exposure_grps, aes(sdrt, n, fill = exposure_category)) +
  geom_col(position = "stack") + base +
  labs(y = "Samples (count)",
       title = "sample-level antibiotic exposures\n(two-day exposure window prior to each fecal sample)") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

# E2c: 100% stacked percentage
e2c <- ggplot(exposure_grps, aes(sdrt, perc, fill = exposure_category)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) + base +
  labs(y = "Samples (%)")

ggsave(file.path(out_dir, "E2b_abx_exposure_count.pdf"), e2b, width = 6.5, height = 4.6)
ggsave(file.path(out_dir, "E2c_abx_exposure_percent.pdf"), e2c, width = 6.5, height = 4.4)
message("wrote results/E2b_abx_exposure_count.pdf, E2c_abx_exposure_percent.pdf")
