# Figure 1 e-h: TaxUMAP of the combined cohort, colour-coded.
#
# The downstream (plotting) half of 162_taxUMAP__Figure1_e_to_h.Rmd. It consumes a
# TaxUMAP embedding that was produced once by the Python tool (see
# taxumap_pipeline_HOWTO.md for how to regenerate taxumap_embedding.csv) and joins
# it to the diet tracker and the diet alpha-diversity to colour the four panels:
#   F1e  day relative to transplant (binned)
#   F1f  daily caloric intake
#   F1g  most-consumed food group (with in-cluster labels)
#   F1h  diet alpha-diversity (Faith PD)
#
# Inputs (all in released_data/): taxumap_embedding.csv, 152_combined_DTB.csv,
# food_group_color_key_final.csv, diet-alpha-diversity.tsv. This is a deterministic
# plotting step, so nothing is cached.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrastr)
  library(RColorBrewer)
  library(viridis)
  library(here)
})

released <- function(f) file.path(Sys.getenv("NUTRITION_DATA", unset = here::here("released_data")), f)
results_dir <- here::here("results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

umap_pt_size <- 0.6
gray_bg <- "gray98"

# A shared theme: clean gray-background scatter with TaxUMAP axis titles, no ticks,
# and a visible legend on the right (the original hid every legend).
theme_taxumap <- function() {
  theme_classic() +
    theme(line = element_blank(),
          legend.position = "right", legend.title = element_text(size = 8),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_text(size = 9, colour = "gray30"),
          panel.background = element_rect(fill = gray_bg),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = gray_bg),
          aspect.ratio = 1)
}
umap_points <- function(mapping) rasterise(geom_point(mapping, alpha = 1, size = umap_pt_size, shape = 16), dpi = 300)

# Inputs ----------------------------------------------------------------------
taxumap   <- read_csv(released("taxumap_embedding.csv"), show_col_types = FALSE)
dtb       <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE)
color_key <- read_csv(released("food_group_color_key_final.csv"), col_types = "ccccc")
faith     <- read_tsv(released("diet-alpha-diversity.tsv"), show_col_types = FALSE) %>% rename(sampleid = 1)

# Per (pid, fdrt): dominant food group (by dehydrated weight) + total calories ---
top_fg <- dtb %>%
  mutate(fgrp1 = str_sub(as.character(Food_code), 1, 1)) %>%
  group_by(pid, fdrt, fgrp1) %>%
  summarise(fg_dehydrated = sum(dehydrated_weight, na.rm = TRUE), .groups = "drop") %>%
  group_by(pid, fdrt) %>%
  slice_max(fg_dehydrated, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(color_key, by = "fgrp1") %>%
  select(pid, fdrt, top_fg_name = shortname, top_fg_color = color)

cal_totals <- dtb %>%
  group_by(pid, fdrt) %>%
  summarise(total_calories = sum(Calories_kcal, na.rm = TRUE), .groups = "drop")

# Combine, attach Faith PD (sampleid = "{pid}d{fdrt}") and TaxUMAP coords. -------
plot_tbl <- top_fg %>%
  left_join(cal_totals, by = c("pid", "fdrt")) %>%
  mutate(sampleid = paste0(pid, "d", fdrt)) %>%
  left_join(faith, by = "sampleid") %>%
  mutate(index_column = paste0("P", pid, "_", fdrt)) %>%
  inner_join(taxumap, by = "index_column")

# F1e: day relative to transplant ---------------------------------------------
# Split the day axis into 3 pre-transplant and 6 post-transplant bins, ordered
# from earliest to latest so the Spectral ramp reads in time order.
bins <- plot_tbl %>% distinct(fdrt) %>%
  mutate(side = fdrt <= 0) %>%
  group_by(side) %>%
  mutate(bin = if (first(side)) as.character(cut_number(fdrt, 3)) else as.character(cut_number(fdrt, 6))) %>%
  ungroup()
bin_levels <- bins %>% arrange(fdrt) %>% pull(bin) %>% unique()

umap_time <- plot_tbl %>%
  left_join(bins %>% select(fdrt, bin), by = "fdrt") %>%
  mutate(bin = factor(bin, levels = bin_levels)) %>%
  arrange(desc(fdrt))

p_e <- ggplot(umap_time) +
  umap_points(aes(taxumap1, taxumap2, colour = bin)) +
  scale_colour_manual(values = brewer.pal(length(bin_levels), "Spectral"), name = "day rel.\nto HCT") +
  labs(x = "TaxUMAP1", y = "TaxUMAP2", title = "Day relative to transplant") +
  theme_taxumap() +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
ggsave(file.path(results_dir, "F1e_taxumap_day.pdf"), p_e, width = 4.2, height = 3.2)

# F1f: caloric intake (sqrt colour spread, axis in 10^3 kcal) ------------------
p_f <- plot_tbl %>%
  mutate(kcal_k = total_calories / 1000) %>%
  ggplot() +
  umap_points(aes(taxumap1, taxumap2, colour = kcal_k)) +
  scale_colour_viridis(option = "plasma", trans = "sqrt", breaks = c(0, 1, 2, 3),
                       name = expression("(x" ~ 10^3 ~ "Kcal)")) +
  labs(x = "TaxUMAP1", y = "TaxUMAP2", title = "Caloric intake") +
  theme_taxumap()
ggsave(file.path(results_dir, "F1f_taxumap_calories.pdf"), p_f, width = 4.2, height = 3.2)

# F1g: most-consumed food group, with in-cluster labels -----------------------
fg_palette <- plot_tbl %>% distinct(top_fg_name, top_fg_color) %>% deframe()
labelled_groups <- c("Meats", "Fruits", "Sweets", "Milk", "Grains")
fg_label_pos <- plot_tbl %>%
  filter(top_fg_name %in% labelled_groups) %>%
  group_by(top_fg_name) %>%
  summarise(taxumap1 = median(taxumap1), taxumap2 = median(taxumap2), .groups = "drop")

p_g <- ggplot(plot_tbl) +
  umap_points(aes(taxumap1, taxumap2, colour = top_fg_name)) +
  # white halo behind, then the group-coloured label on top, so it stays legible
  # against its own cluster.
  geom_text(data = fg_label_pos, aes(taxumap1, taxumap2, label = top_fg_name),
            colour = "white", fontface = "bold", size = 4.4, show.legend = FALSE) +
  geom_text(data = fg_label_pos, aes(taxumap1, taxumap2, label = top_fg_name, colour = top_fg_name),
            fontface = "bold", size = 4, show.legend = FALSE) +
  scale_colour_manual(values = fg_palette, name = NULL) +
  labs(x = "TaxUMAP1", y = "TaxUMAP2", title = "Most consumed food group") +
  theme_taxumap() +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
ggsave(file.path(results_dir, "F1g_taxumap_foodgroup.pdf"), p_g, width = 4.6, height = 3.2)

# F1h: diet alpha-diversity (Faith PD, axis in 10^3) --------------------------
p_h <- plot_tbl %>%
  filter(!is.na(faith_pd)) %>%
  mutate(faith_k = faith_pd / 1000) %>%
  ggplot() +
  umap_points(aes(taxumap1, taxumap2, colour = faith_k)) +
  scale_colour_viridis(breaks = c(0, 1, 2), name = expression("(x" ~ 10^3 ~ ")")) +
  labs(x = "TaxUMAP1", y = "TaxUMAP2", title = expression("Diet" ~ alpha * "-diversity")) +
  theme_taxumap()
ggsave(file.path(results_dir, "F1h_taxumap_diversity.pdf"), p_h, width = 4.2, height = 3.2)

message("TaxUMAP panels F1e-h written to results/.")
