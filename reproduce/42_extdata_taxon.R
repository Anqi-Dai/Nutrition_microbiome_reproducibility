# Extended Fig. E7 (human taxon family): E7a, E7b and E7e.
#
# Refactor of R63_Enterococcus_asv1_outcome.Rmd (E7a, E7e) and
# 178_new_F4__code_for_Figure_4.Rmd (E7b).
#   E7a  food-group effect-size heatmap across the prevalent genera, rows
#        clustered (dendrogram) by their posterior fraction-positive profile
#   E7b  genus-abundance vs alpha-diversity Spearman correlations across the same
#        genera (the full-set companion to F4a)
#   E7e  marginal E. faecium (asv_1) CLR over each food group, split by
#        antibiotic exposure, from the cached asv_1 fit's conditional effects
#
# Plot-only for E7a (reads the per-genus summary cached by 40); E7b reuses the
# shared Spearman helper; E7e reloads the cached asv_1 brms fit for its
# conditional effects.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages({
  library(dendextend)
  library(ggdendro)
})

key <- food_key()
replacement_dictionary <- setNames(key$shortname, key$fg1_name)

# Panel proportions measured off the assembled Extended Fig. 7 artwork: a 142 x 217 pt
# heatmap and a 40 x 211 pt bar panel, both on a 6.6 pt genus row. `hairline` (0.5 pt)
# and `genus_text_pt` (7 pt) come from the helper, shared with F4a in 41. The rest of
# the axis text is set to the genus size so each panel reads at one size.
heatmap_ratio <- 217 / 142
spearman_ratio <- 211 / 40

# E7a: genus heatmap + dendrogram -----------------------------------------------
level_order <- rev(c(
  "abx * Sweets", "Sweets",
  "abx * Grains", "Grains",
  "abx * Milk", "Milk",
  "abx * Eggs", "Eggs",
  "abx * Legumes", "Legumes",
  "abx * Meats", "Meats",
  "abx * Fruits", "Fruits",
  "abx * Oils", "Oils",
  "abx * Vegetables", "Vegetables"))

# Keep only the food-group terms; assign tiered significance marks from the
# nested credible intervals (*** 99% / ** 97% / * 94%), checking the strictest
# threshold first and cascading down.
post_summary <- read_csv(cache_path("R63_genus_clr_all_models_results.csv"),
                         show_col_types = FALSE) |>
  filter(str_detect(term, "fg_")) |>
  mutate(
    clean_term = term |>
      str_replace("empiricalTRUE$", "abx") |>
      str_remove_all("empiricalFALSE:|avg_intake_|TRUE$") |>
      str_replace_all(replacement_dictionary) |>
      str_replace("empiricalTRUE:", "abx * ") |>
      str_replace_all("_", " "),
    mark = case_when(
      Q0.5 > 0 ~ "***",
      Q99.5 < 0 ~ "***",
      Q1.5 > 0 ~ "**",
      Q98.5 < 0 ~ "**",
      Q3 > 0 ~ "*",
      Q97 < 0 ~ "*",
      .default = ""),
    clean_term = factor(clean_term, levels = level_order))

# Cluster genera (rows) on their fraction-positive profile across the food terms.
post_matrix <- post_summary |>
  filter(!str_detect(term, "intensity")) |>
  select(outcome, term, fraction_positive) |>
  spread("term", "fraction_positive") |>
  column_to_rownames("outcome")


plot_hc_w_node_labels <- function(hc, ...){
  plot(hc, ...)
  dend <- as.dendrogram(hc)
  xy <- get_nodes_xy(dend)
  heights <- get_nodes_attr(dend, "height")
  internal <- heights > 0
  text(
    x = xy[internal, 1],
    y = xy[internal, 2],
    labels = seq_len(sum(internal)),
    pos = 3,
    cex = 0.7,
    col = "red"
  )
}
adjust_aesthetics = FALSE
genus_distance <- dist(post_matrix)
hc <- hc_raw <- hclust(genus_distance, method = "complete")
if (adjust_aesthetics) {
  labels(hc) <- paste0(1:length(hc$labels), " " , labels(hc))
  plot_hc_w_node_labels(hc_raw)
}

# Row order is rotated towards the printed panel. Rotation only chooses which child of
# each split is drawn first, which is arbitrary in a dendrogram, so every cluster and
# every merge height is left untouched and only the reading order changes.
#
# The printed row order, top to bottom. Genera the printed panel does not carry take
# their nearest neighbour's target position, so they land beside the genus they cluster
# with instead of dragging the sort around.
printed_row_order <- c(
  "Streptococcus", "Lactococcus", "Faecalimonas", "Romboutsia", "Blautia",
  "Lachnoclostridium", "Flavonifractor", "Erysipelatoclostridium", "Anaerostipes",
  "Eggerthella", "Ruminococcus", "Escherichia", "Fusicatenibacter", "Bifidobacterium",
  "Akkermansia", "Parabacteroides", "Bacteroides", "Dorea", "Clostridium",
  "Anaerobutyricum", "Schaalia", "Actinomyces", "Scardovia", "Intestinibacter",
  "Staphylococcus", "Rothia", "Granulicatella", "Veillonella", "Atopobium",
  "Lactobacillus", "Enterococcus")

# Target positions run bottom to top, since the first level of the y factor is the
# bottom row of the panel.
target_position <- setNames(as.numeric(seq_along(printed_row_order)), rev(printed_row_order))
distance_matrix <- as.matrix(genus_distance)
for (genus in setdiff(hc$labels, names(target_position))) {
  nearest <- names(sort(distance_matrix[genus, names(target_position)]))[1]
  target_position[genus] <- target_position[[nearest]] + 0.5
}
target_position <- rank(target_position[hc$labels])

# Pick the best of the 2^(n-1) rotations exactly, rather than hand-tuning branch flips:
# a dendrogram rotation is a free choice of child order at each split, so the arrangement
# that minimises total row displacement from the printed order can be found by working up
# from the leaves. For a subtree placed at a known starting row the two child orders are
# independent, so each is solved once and cached.
best_rotation <- local({
  memo <- new.env(hash = TRUE)
  subtree_size <- function(node) if (node < 0) 1L else
    subtree_size(hc$merge[node, 1]) + subtree_size(hc$merge[node, 2])
  solve_node <- function(node, start) {
    key <- paste(node, start)
    cached <- memo[[key]]
    if (!is.null(cached)) return(cached)
    result <- if (node < 0) {
      genus <- hc$labels[-node]
      list(cost = abs(start - target_position[[genus]]), order = genus)
    } else {
      left <- hc$merge[node, 1]; right <- hc$merge[node, 2]
      as_is <- list(solve_node(left, start), solve_node(right, start + subtree_size(left)))
      flipped <- list(solve_node(right, start), solve_node(left, start + subtree_size(right)))
      pick <- if (sum(sapply(as_is, `[[`, "cost")) <= sum(sapply(flipped, `[[`, "cost"))) as_is else flipped
      list(cost = sum(sapply(pick, `[[`, "cost")),
           order = c(pick[[1]]$order, pick[[2]]$order))
    }
    assign(key, result, envir = memo)
    result
  }
  solve_node(nrow(hc$merge), 1)
})

hc <- rotate(hc, best_rotation$order)
dendro_data <- dendro_data(hc, type = "rectangle")
ordered_genera <- hc$labels[hc$order]

# ordered_genera[1] is the bottom row. Enterococcus is asked for the bottom and
# Streptococcus/Lactococcus for the top two; the optimisation puts them there on this
# tree, but a future refit could move them into one cluster, where no rotation can
# reach both ends. Stop rather than print a panel that quietly drops the request.
stopifnot(identical(ordered_genera, best_rotation$order),
          ordered_genera[1] == "Enterococcus",
          tail(ordered_genera, 2) == c("Lactococcus", "Streptococcus"))

dendro_plot <- ggplot(segment(dendro_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0.5)) +
  scale_x_continuous(expand = c(0, 0.5)) +
  theme_dendro()

heatmap_data <- post_summary |>
  mutate(outcome = factor(outcome, levels = ordered_genera))

heatmap_plot <- heatmap_data |>
  # Blank the non-significant tiles so the credible effects read clearly.
  mutate(Estimate = if_else(!str_detect(mark, "\\*"), 0, Estimate)) |>
  ggplot(aes(x = clean_term, y = outcome, fill = Estimate)) +
  geom_tile(color = "black", linewidth = hairline, width = 0.95, height = 0.95) +
  geom_text(aes(label = mark), nudge_y = -0.15, size = 1.5, color = "black") +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "firebrick", midpoint = 0) +
  labs(fill = "Effect size\n(posterior mean)") +
  theme_minimal(base_size = 11) +
  theme(
    aspect.ratio = heatmap_ratio,
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown(face = "italic", size = genus_text_pt),
    axis.text.x = element_markdown(angle = 45, hjust = 1, size = genus_text_pt),
    axis.ticks = element_line(linewidth = hairline),
    legend.position = "left",
    legend.title = element_text(size = genus_text_pt),
    legend.text = element_text(size = genus_text_pt),
    panel.background = element_rect(colour = "white", fill = "white"),
    panel.grid = element_blank())

final_figure <- plot_grid(heatmap_plot, dendro_plot, ncol = 2,
                          rel_widths = c(1, 0.1), align = "h", axis = "tb")
save_panel(final_figure, "E7a_genus_heatmap.pdf", width = 118, height = 99)
message("E7a done")

# E7b: Spearman correlations across the prevalent genera ------------------------
# Same genus set as the heatmap above and as F4a, from the shared prevalence filter
# in the helper. Keep the FDR-significant ones, ordered by signed rho. The set comes
# from the filter rather than from the cached 40 results, so E7b does not depend on
# the per-genus models having been fitted.
res <- genus_diversity_spearman()

selected_bars <- res %>%
  filter(sig05 == "FDR < 0.05") %>%
  arrange(desc(rho))

correbar_all <- selected_bars %>%
  mutate(genus = factor(genus, levels = selected_bars$genus)) %>%
  ggplot(aes(x = genus, y = 0, xend = genus, yend = rho, color = Correlation)) +
  # bar thickness in points, so it stays proportionate to the 6.6 pt genus row
  geom_segment(linewidth = 4.5 / .pt) +
  # Axis titled with the symbol alone and only three breaks: at 40 pt wide the panel
  # has no room for five tick labels or a spelled-out title, and this is what the
  # assembled figure shows (the "Spearman correlation with alpha-diversity" wording
  # sits beside the panel as artwork, not as an axis title).
  labs(x = "", y = expression(rho)) +
  scale_y_continuous(breaks = c(-0.2, 0, 0.4), labels = c("-0.2", "0", "0.4")) +
  scale_color_jco() +
  coord_flip() +
  theme_classic(base_size = 10) +
  theme(aspect.ratio = spearman_ratio,
        axis.text = element_text(size = genus_text_pt),
        # the panel is only 40 pt wide, so the three tick labels are set a step
        # below the genus names to keep them from touching, as in the artwork
        axis.text.x = element_text(size = genus_text_pt - 1),
        axis.title = element_text(size = genus_text_pt),
        legend.position = "none",
        axis.line = element_line(linewidth = hairline),
        axis.ticks = element_line(linewidth = hairline),
        axis.text.y = element_markdown(face = "italic", size = genus_text_pt))

save_panel(correbar_all, "E7b_genus_diversity_spearman.pdf", width = 52, height = 83)
message("E7b done")

# E7e: marginal E. faecium CLR per food group, by antibiotic exposure -----------
asv1_fit <- readRDS(cache_path("R63_fit_asv1.rds"))

# Conditional effects for the food-group x antibiotics interactions. Stretch the
# per-100g model units back to grams for the axes.
raw_conditional <- conditional_effects(asv1_fit, surface = TRUE)
condi_dat <- raw_conditional[str_detect(names(raw_conditional), ":empirical")] |>
  bind_rows(.id = "grp") |>
  mutate(grp = str_replace(grp, "\\:empirical", "")) |>
  mutate(effect1__ = effect1__ * 100,
         across(starts_with("fg_"), ~ .x * 100)) |>
  left_join(key |> select(grp = fg1_name, shortname), by = "grp") |>
  mutate(shortname = factor(shortname, levels = rev(c(
    "Vegetables", "Oils", "Fruits", "Meats",
    "Legumes", "Eggs", "Milk", "Grains", "Sweets"))))

conditional_plots <- condi_dat |>
  ggplot() +
  geom_smooth(aes(x = effect1__, y = estimate__, ymin = lower__, ymax = upper__,
                  fill = effect2__, color = effect2__),
              stat = "identity", alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual("broad-spectrum antibiotics", values = abx_palette,
                    labels = c("not exposed", "exposed")) +
  scale_colour_manual("broad-spectrum antibiotics", values = abx_palette,
                      labels = c("not exposed", "exposed")) +
  facet_wrap(~ shortname, nrow = 3, scales = "free_x") +
  labs(y = "Predicted CLR(E. faecium)", x = "Food group consumed (grams)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8, face = "bold"),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.key = element_rect(fill = alpha("white", 0)),
        aspect.ratio = 1,
        strip.background = element_rect(color = "white", fill = "#ffecdc", linewidth = 1.5, linetype = "solid")) +
  guides(fill = guide_legend(direction = "vertical"), color = guide_legend(direction = "vertical"))

save_panel(conditional_plots, "E7e_efaecium_marginal.pdf", width = 150, height = 160)
message("E7e done")
