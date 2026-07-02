# Extended Data E6i: cumulative incidence of hospital discharge following neutrophil
# engraftment, by diet-pattern cluster (RESTRICTED). Companion to R09/R10's clinical
# outcome analysis; the discharge outcome is landmarked at engraftment.
#
# Time = tLOS = tdischarge - tengraftment (days after neutrophil engraftment), event =
# LOS (1 if discharged). The cohort is the patients with a valid engraftment landmark
# (tengraftment present); every one of them is eventually discharged, so the cumulative
# incidence is 1 - KM and rises to 1. Cluster comes from the cleaned restricted
# df_main_clinical_outcome.rds, so this is a restricted panel that skips when absent.
#
# Two pieces:
#   - curves: cumulative incidence of discharge by cluster (Cluster 1 navy, Cluster 2
#     gold), with the published at-risk / cumulative-event table below.
#   - adjusted HR: a Cox model of discharge on cluster adjusted for the duration of
#     broad-spectrum antibiotic exposure (day_exposed, day -7..12), graft source,
#     conditioning intensity and GVHD prophylaxis. Cluster 2 vs Cluster 1 HR annotated.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({
  library(survival)
  library(ggsurvfit)
  library(patchwork)
})

df_file <- "df_main_clinical_outcome.rds"
if (!has_restricted(df_file)) {
  message("E6i skipped: restricted df_main not found (", restricted(df_file), ").")
  message("This panel needs the diet-pattern cluster and discharge landmark; place the ",
          "cleaned df_main_clinical_outcome.rds in restricted_data/.")
  quit(save = "no", status = 0)
}

results_dir <- here::here("results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cluster_cols <- c("Cluster 1" = "darkslateblue", "Cluster 2" = "darkgoldenrod2")

# Cohort: patients with a valid engraftment landmark (Cluster 1 = 111, Cluster 2 = 57).
# Cluster 1 is the reference so the reported HR is Cluster 2 vs Cluster 1.
df <- read_rds(restricted(df_file)) |>
  filter(!is.na(tengraftment)) |>
  mutate(modal_diet = factor(modal_diet, levels = c("Cluster 1", "Cluster 2")),
         source     = factor(source),
         intensity  = factor(intensity),
         gvhd_ppx   = factor(gvhd_ppx))

message("E6i cohort: ", nrow(df), " patients (",
        paste(names(table(df$modal_diet)), table(df$modal_diet), sep = " = ", collapse = ", "), ").")

# Adjusted Cox model: discharge hazard by cluster, adjusted for broad-spectrum abx
# exposure duration, graft source, conditioning intensity and GVHD prophylaxis.
fit_cox <- coxph(Surv(tLOS, LOS) ~ modal_diet + day_exposed + source + intensity + gvhd_ppx,
                 data = df)
cox_row <- summary(fit_cox)$coefficients["modal_dietCluster 2", ]
hr <- cox_row["exp(coef)"]
pv <- cox_row["Pr(>|z|)"]
hr_lab <- sprintf("adjusted HR=%.2f\np=%.3f", hr, pv)
message("E6i adjusted HR (Cluster 2 vs 1) = ", sprintf("%.2f", hr),
        ", p = ", sprintf("%.3f", pv))

# Cumulative incidence of discharge (type = "risk" -> 1 - KM). The curve reaches 1
# because every landmark patient is eventually discharged (no competing event).
xbreaks <- seq(0, 40, 10)
curve <- survfit2(Surv(tLOS, LOS) ~ modal_diet, data = df) |>
  ggsurvfit(type = "risk", linewidth = 0.9) +
  scale_colour_manual(values = cluster_cols, labels = names(cluster_cols)) +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0)) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 1.02)) +
  annotate("text", x = 21, y = 0.72, hjust = 0, label = hr_lab, size = 4.2) +
  labs(x = "days after neutrophil engraftment",
       y = "hospital discharge\ncumulative incidence") +
  theme(legend.position = c(0.46, 0.26))

# Risk table below the curve, laid out to match the published panel: one block per
# cluster (header + At Risk + cumulative Events rows), coloured to its curve. Built
# as its own ggplot so the text can be cluster-coloured (ggsurvfit's own risktable
# only colours when grouped by statistic, not by strata). survfit's per-time n.risk
# and n.event at the tick days give the two rows; events are cumulated.
sf  <- survfit(Surv(tLOS, LOS) ~ modal_diet, data = df)
smy <- summary(sf, times = xbreaks, extend = TRUE)
tab <- tibble(cluster = sub("modal_diet=", "", as.character(smy$strata)),
              time = smy$time, n.risk = smy$n.risk, event = smy$n.event) |>
  group_by(cluster) |> mutate(event = cumsum(event)) |> ungroup()

# Fixed row heights (top to bottom): header, At Risk, Events, per cluster.
rows <- tribble(~cluster,     ~stat,      ~y,
                "Cluster 1",  "At Risk",  4.8,
                "Cluster 1",  "Events",   4.0,
                "Cluster 2",  "At Risk",  1.8,
                "Cluster 2",  "Events",   1.0)
hdr  <- tibble(cluster = c("Cluster 1", "Cluster 2"), y = c(5.7, 2.7))
vals <- tab |>
  pivot_longer(c(n.risk, event), names_to = "stat", values_to = "val") |>
  mutate(stat = if_else(stat == "n.risk", "At Risk", "Events")) |>
  left_join(rows, by = c("cluster", "stat"))

# Same x window as the curve so the values sit under the axis ticks; the row/cluster
# labels sit at negative x (in the shared left margin) with clipping off.
risktable <- ggplot() +
  geom_text(data = hdr,  aes(-8,   y, label = paste0(cluster, ":"), colour = cluster),
            hjust = 0, fontface = "bold", size = 3.8) +
  geom_text(data = rows, aes(-1.5, y, label = stat, colour = cluster), hjust = 1, size = 3.5) +
  geom_text(data = vals, aes(time, y, label = val,  colour = cluster), size = 3.7) +
  scale_colour_manual(values = cluster_cols, guide = "none") +
  scale_x_continuous(breaks = xbreaks) +
  scale_y_continuous(limits = c(0.4, 6.2)) +
  coord_cartesian(xlim = c(0, 40), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 10, 2, 8))

panel <- curve / risktable + plot_layout(heights = c(1, 0.46))
ggsave(file.path(results_dir, "E6i_discharge_cuminc.pdf"), panel, width = 6, height = 6)

message("Wrote results/E6i_discharge_cuminc.pdf")
