# Data S6: per-patient timecourse of daily calories, diet Faith diversity, and
# fecal alpha-diversity (RESTRICTED).
#
# Ported from 085_each_pt_timecourse__code_for_Figure_S2.Rmd. One multi-page A4
# PDF: an enlarged example patient sits top-left on page 1, every other patient is
# a small panel, each page carries one shared "Transplant day" x-axis label.
#
# Restricted because it needs the per-patient engraftment day
# (085_engraftment_day_annot.csv, the green dashed line): PHI-free but not cleared
# for public release. The other three inputs (DTB, META, diet Faith PD) are all in
# released_data. Skips cleanly when restricted_data/ is absent.

source(here::here("reproduce", "human", "_human_helpers.R"))

engraft_file <- "085_engraftment_day_annot.csv"
if (!has_restricted(engraft_file)) {
  message("Data S6 skipped: restricted input not found (", restricted(engraft_file), ").")
  message("This figure needs the per-patient engraftment day, which is not ",
          "publicly released. Point RESTRICTED_DATA at it or place it in restricted_data/.")
  quit(save = "no", status = 0)
}

cfg <- list(
  example_pid     = "P1",         # shown enlarged, with the detailed Diet:/Stool: labels
  scale_factor    = 92,           # 4600/50 -- brings stool diversity onto the calorie axis
  y_limit         = 4600,         # shared upper bound so every panel is comparable
  pt_line_size    = 0.2,
  axis_line_thick = 0.2,
  diet_color      = "#E41A1C",
  stool_color     = "blue",
  faith_color     = "black",
  panels_per_page = 40,           # 5 columns x 8 rows on an A4 sheet
  ncol            = 5
)

out_file <- here::here("results", "dataS6_all_patients_timecourse.pdf")

# Load and reshape. Each dataset is split into a per-patient list so a panel is
# built by looking up one patient's name.
dtb  <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE)
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) |> split(~ pid)

# daily calories: sum every food record within a patient-day
day_calori <- dtb |>
  group_by(pid, fdrt) |>
  summarise(daycal = sum(Calories_kcal), .groups = "drop") |>
  split(~ pid)

# diet Faith diversity: the released diet-alpha-diversity.tsv packs patient and
# day into the first (unnamed) column as "<pid>d<day>"
faith_path <- released("diet-alpha-diversity.tsv")
faith <- read_tsv(faith_path, show_col_types = FALSE) |>
  rename(sampleid = 1) |>
  separate(sampleid, into = c("pid", "fdrt"), sep = "d") |>
  mutate(fdrt = as.numeric(fdrt)) |>
  split(~ pid)

# patients carried through the figure are those that actually have stool data
pids <- names(meta)

# per-patient annotation: engraftment day (restricted), diet-day and stool-sample counts
engraftment <- read_csv(restricted(engraft_file), show_col_types = FALSE) |>
  filter(pid %in% pids)

diet_days <- dtb |>
  filter(pid %in% pids) |>
  distinct(pid, fdrt) |>
  count(pid, name = "d_days")

stool_days <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE) |>
  filter(pid %in% pids) |>
  distinct(pid, sdrt) |>
  count(pid, name = "s_days")

annot <- engraftment |>
  full_join(diet_days, by = "pid") |>
  full_join(stool_days, by = "pid") |>
  split(~ pid)

message("Faith pd range:  ", paste(range(read_tsv(faith_path, show_col_types = FALSE)$faith_pd), collapse = " - "))
message("Daily cal range: ", paste(range(day_calori |> bind_rows() |> pull(daycal)), collapse = " - "))

# One function builds any patient's panel. detailed = TRUE gives the enlarged
# example look (axis titles, "Diet:"/"Stool:" prefixes); FALSE the bare grid panel.
make_panel <- function(pid, detailed = FALSE) {

  diet_label  <- if (detailed) "Diet: {d_days} days"     else "{d_days} days"
  stool_label <- if (detailed) "Stool: {s_days} samples" else "{s_days} samples"

  p <- ggplot() +
    geom_line(data = day_calori[[pid]], aes(fdrt, daycal),
              linewidth = cfg$pt_line_size, colour = cfg$diet_color) +
    geom_line(data = faith[[pid]], aes(fdrt, faith_pd),
              linewidth = cfg$pt_line_size, colour = cfg$faith_color) +
    geom_line(data = meta[[pid]], aes(sdrt, simpson_reciprocal * cfg$scale_factor),
              linewidth = cfg$pt_line_size, colour = cfg$stool_color) +
    geom_vline(data = annot[[pid]], aes(xintercept = engraftment_day),
               colour = "forestgreen", linewidth = cfg$axis_line_thick, linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "darkgray", linewidth = cfg$axis_line_thick) +
    geom_text(data = annot[[pid]], inherit.aes = FALSE,
              aes(Inf, Inf, label = str_glue(diet_label)),
              hjust = 1.6, vjust = 1, size = 2, colour = cfg$diet_color) +
    geom_text(data = annot[[pid]], inherit.aes = FALSE,
              aes(Inf, Inf, label = str_glue(stool_label)),
              hjust = 1.3, vjust = 2.5, size = 2, colour = cfg$stool_color) +
    theme_pubr()

  if (detailed) {
    p +
      scale_y_continuous(
        name = "Daily calories, Faith div", limits = c(0, cfg$y_limit),
        sec.axis = sec_axis(~ . / cfg$scale_factor,
                            name = expression(Fecal ~ alpha ~ diversity))) +
      labs(x = "Transplant day", y = "", title = pid) +
      theme(axis.text          = element_text(size = 5),
            axis.title         = element_text(size = 8),
            plot.title         = element_text(size = 8),
            axis.title.y       = element_text(colour = cfg$diet_color),
            axis.title.y.right = element_text(colour = cfg$stool_color),
            axis.line          = element_line(colour = "black", linewidth = cfg$axis_line_thick),
            axis.ticks         = element_line(colour = "black", linewidth = cfg$axis_line_thick),
            aspect.ratio       = 1/2)
  } else {
    # small grid panels drop axis titles and the left y numbers to cut clutter;
    # tick marks stay so each panel keeps its own scale
    p +
      scale_y_continuous(name = "", limits = c(0, cfg$y_limit),
                         sec.axis = sec_axis(~ . / cfg$scale_factor)) +
      labs(x = "Transplant day", y = "", title = pid) +
      theme(axis.text          = element_text(size = 7),
            axis.title         = element_text(size = 5),
            plot.title         = element_text(size = 8),
            axis.title.x       = element_blank(),
            axis.text.y        = element_blank(),
            axis.title.y       = element_text(colour = cfg$diet_color),
            axis.title.y.right = element_text(colour = cfg$stool_color),
            axis.line          = element_line(colour = "black", linewidth = cfg$axis_line_thick),
            axis.ticks         = element_line(colour = "black", linewidth = cfg$axis_line_thick),
            aspect.ratio       = 1/2)
  }
}

small_panels  <- pids |> set_names() |> map(make_panel, detailed = FALSE)
example_panel <- make_panel(cfg$example_pid, detailed = TRUE)

# Assembles one A4 page: bold title, panel grid, one shared bottom x label.
# Page 1 reserves a top-left 2x2 block for the enlarged example.
xlab_strip <- ggdraw() + draw_label("Transplant day", size = 8)

make_page <- function(panels, page_title, big_panel = NULL) {

  title <- ggdraw() +
    draw_label(page_title, fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))

  if (is.null(big_panel)) {
    body <- plot_grid(plotlist = panels, ncol = cfg$ncol)
  } else {
    # the example fills a top-left 2x2; the first 6 small panels sit to its right,
    # the remaining 30 fill a normal grid below
    top_right <- plot_grid(plotlist = panels[1:6], ncol = 3)
    top_band  <- plot_grid(big_panel, top_right, ncol = 2, rel_widths = c(2, 3))
    bottom    <- plot_grid(plotlist = panels[7:36], ncol = cfg$ncol)
    body      <- plot_grid(top_band, bottom, ncol = 1, rel_heights = c(2, 5))
  }

  plot_grid(title, body, xlab_strip, ncol = 1, rel_heights = c(0.08, 1, 0.05)) +
    theme(plot.margin = unit(c(1, 1, 2, 1), "cm"))
}

# Page 1 holds 36 small panels because the example eats 4 cells; later pages hold
# a full 40. Chunking is computed so adding or losing patients needs no re-editing.
n            <- length(pids)
page1_n      <- cfg$panels_per_page - 4
later_idx    <- seq(page1_n + 1, n)
later_chunks <- split(later_idx, ceiling(seq_along(later_idx) / cfg$panels_per_page))

pages <- c(
  list(make_page(small_panels[1:page1_n], "Data S6", big_panel = example_panel)),
  map(later_chunks, ~ make_page(small_panels[.x], "Data S6 (Continued)"))
)

if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)
pdf(out_file, width = 210 / 25.4, height = 297 / 25.4)  # A4, mm -> inches
walk(pages, print)
dev.off()

message("Wrote ", length(pages), "-page figure to ", out_file)
