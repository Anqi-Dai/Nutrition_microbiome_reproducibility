# E2a: per-patient daily antibiotic-exposure heatmap (RESTRICTED).
#
# Ported from R21_abx_groups__code_for_Figure_S4.Rmd. Unlike E2b/E2c (26), which
# only need the two-day-prior-to-sample exposure window shipped in Data_S4, this
# panel needs each patient's full daily antibiotic time course. That table
# (R21_meds_updated_all_medication_classified.csv: one row per drug course with
# start/stop transplant days and a per-study drug category) is PHI-free but not
# cleared for public release, so it lives in the gitignored restricted_data/ tier.
# The script skips cleanly when that folder is absent.
#
# The stool-sample day markers (thick black boxes) come from the released META.

source(here::here("reproduce", "human", "_human_helpers.R"))

meds_file <- "R21_meds_updated_all_medication_classified.csv"
if (!has_restricted(meds_file)) {
  message("E2a skipped: restricted input not found (", restricted(meds_file), ").")
  message("This panel needs the full daily antibiotic time course, which is not ",
          "publicly released. Point RESTRICTED_DATA at it or place it in restricted_data/.")
  quit(save = "no", status = 0)
}

# One row per drug course: pid, drug_name_clean, drug_route ({drug}__{route}),
# start/stop transplant day, and the per-study category. drug_route is carried
# only to keep the course grouping unique (it was `together` in R21).
meds_updated <- read_csv(restricted(meds_file), show_col_types = FALSE)

# Expand each course to one row per exposed day, then collapse to one exposure
# category per patient-day. Broad-spectrum and fluoroquinolone on the same day ->
# "both"; otherwise the strongest single class present that day.
meds_expanded <- meds_updated |>
  distinct() |>
  group_by(pid, drug_name_clean, drug_route, startday, stopday, drug_category_for_this_study) |>
  reframe(day = seq(startday, stopday)) |>
  group_by(pid, day) |>
  summarise(
    exposure_category = case_when(
      any(drug_category_for_this_study == "broad_spectrum") &
        any(drug_category_for_this_study == "fluoroquinolones") ~ "both",
      any(drug_category_for_this_study == "broad_spectrum")     ~ "broad_spectrum",
      any(drug_category_for_this_study == "fluoroquinolones")   ~ "fluoroquinolones",
      any(drug_category_for_this_study == "other_antibacterials") ~ "other_antibacterials",
      TRUE ~ "not_antibacterial"
    ), .groups = "drop"
  )

# Released META supplies which patient-days have a stool sample.
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) |>
  select(pid, sdrt)

df <- meds_expanded

# Order patients by the first broad-spectrum day (patients with none sink to the
# bottom, keyed at day 40).
first_broad_spectrum <- df |>
  filter(exposure_category == "broad_spectrum") |>
  group_by(pid) |>
  summarize(FirstBroadSpectrumDay = min(day), .groups = "drop") |>
  full_join(df |> distinct(pid), by = "pid") |>
  mutate(FirstBroadSpectrumDay = if_else(is.na(FirstBroadSpectrumDay), 40, FirstBroadSpectrumDay))

sorted_patient_ids <- first_broad_spectrum |>
  arrange(FirstBroadSpectrumDay) |>
  pull(pid)

df$pid <- factor(df$pid)

# Fill in the non-antibacterial days so the grid is complete, flag stool-sample
# days, and lock the patient order.
df_complete <- df |>
  complete(pid, day = full_seq(day, 1), fill = list(exposure_category = "not_antibacterial")) |>
  mutate(has_stool = paste(pid, day) %in% paste(meta$pid, meta$sdrt)) |>
  mutate(pid = factor(pid, levels = sorted_patient_ids))

color_palette <- c(
  "broad_spectrum"       = "#ef7f7f",
  "fluoroquinolones"     = "#b3bce2",
  "both"                 = "#6a4c95",
  "not_antibacterial"    = "white",
  "other_antibacterials" = "#f5c5e1"
)

heatmap <- ggplot(df_complete, aes(x = day, y = pid)) +
  geom_tile(aes(fill = exposure_category),
            color = "black", linewidth = 0.01, show.legend = FALSE) +
  geom_tile(data = df_complete |> filter(has_stool),
            aes(color = has_stool), fill = NA, linewidth = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = color_palette, na.value = NA) +
  scale_color_manual(values = c("TRUE" = "black")) +
  scale_x_continuous(breaks = c(0, 20, 40)) +
  labs(title = "", x = "Transplant day", y = "Patient ID", fill = "Exposure Category") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x  = element_text(size = 11, angle = 0, hjust = 1, vjust = 0.5),
    axis.text.y  = element_text(size = 6,  angle = 0, hjust = 1, vjust = 0.5),
    axis.title   = element_text(size = 11),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 8)
  )

save_panel(heatmap, "E2a_abx_exposure_heatmap.pdf", width = 100, height = 250)

message("E2a exposure heatmap written to results/ (",
        nlevels(df_complete$pid), " patients).")
