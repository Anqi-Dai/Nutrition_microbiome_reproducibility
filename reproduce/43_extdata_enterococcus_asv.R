# Extended Fig. E7 c,d: Enterococcus ASV / species summary, cleaned up from
# 03-summary-figure.Rmd. Two panels:
#
#   E7c  E7c_asv_distribution.pdf -- stacked barplot of 16S Enterococcus ASV
#        composition across the fecal samples whose Enterococcus relative abundance
#        is >= 0.5% (named top-10 ASVs, "Other Enterococcus", "not Enterococcus").
#   E7d  E7d_species_table.pdf -- per-ASV NCBI RefSeq species matches with their
#        metagenomic (MGX) prevalence (n / 333) and abundance distribution.
#
# Inputs are plain CSVs in released_data/ (no phyloseq):
#   171_16S_enterococcus_asv_relab.csv   sampleid, asv, species, relab  (16S
#                                        Enterococcus ASVs, fraction of the whole
#                                        community; the 16S species is carried per
#                                        ASV, so species-level composition -- used
#                                        for the E7c sample order -- is just an
#                                        aggregation of this one table)
#   mgx_enterococcus_species_relab.csv   sample, species, relab (fraction of the
#                                        whole MGX community; Enterococcus only)
# These were exported once from the original phyloseq objects.
#
# The ASV -> species cluster grouping (Enterococcus species sharing an identical
# V4V5 amplicon with each ASV) comes from the sytycs amplicon analysis, which is not
# shipped; it is a fixed 4-ASV lookup and is hard-coded below. The MGX prevalences
# it annotates are recomputed here and reproduce the published panel exactly
# (faecium 171, faecalis 114, gallinarum 33, casseliflavus 31, ...).

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

released <- function(f) here::here("released_data", f)
out <- function(f) here::here("results", f)
if (!dir.exists(here::here("results"))) dir.create(here::here("results"))

N_MGX <- 333L   # total MGX samples (prevalence denominator; fixed released cohort)

# ASV colours (the Enterococcus ASV palette; was colors.R, inlined here)
asv_cols <- c(asv_1 = "#46ACC8", asv_33 = "darkgreen", asv_67 = "brown4",
              asv_19 = "lightgreen", asv_12 = "gold", asv_33138 = "purple4",
              asv_199390 = "red", asv_233 = "#235368", asv_3179 = "#D39C6D",
              asv_16631 = "pink3", asv_211112 = "navy", asv_16968 = "lightgrey")
# sample ordering for E7c. "olo" (optimal leaf ordering on a Bray-Curtis distance,
# via the seriation package) is microViz comp_barplot's own default ordering and
# best matches the published panel; it needs the lightweight `seriation` package and
# falls back to plain Bray hclust if that is not installed. Alternatives: "hclust",
# "pcoa1" (smooth monotonic gradient), "entdesc", "asv1".
ORDER <- Sys.getenv("BAR_ORDER", "olo")
if (ORDER == "olo" && !requireNamespace("seriation", quietly = TRUE)) {
  message("seriation not installed; falling back to hclust ordering")
  ORDER <- "hclust"
}

asv_pal <- c(asv_cols,
             "Other Enterococcus" = "gray60",
             "not Enterococcus"   = "#FCF8D4")   # pale yellow

# =====================================================================
# E7c: stacked barplot of 16S Enterococcus ASV composition
# =====================================================================
e16 <- read_csv(released("171_16S_enterococcus_asv_relab.csv"), show_col_types = FALSE)

top10 <- e16 |> group_by(asv) |> summarise(tot = sum(relab), .groups = "drop") |>
  slice_max(tot, n = 10) |> pull(asv)

ent_long <- e16 |>
  mutate(taxlev = ifelse(asv %in% top10, asv, "Other Enterococcus")) |>
  group_by(sampleid, taxlev) |> summarise(relab = sum(relab), .groups = "drop")

sample_ent <- ent_long |> group_by(sampleid) |>
  summarise(ent_total = sum(relab), .groups = "drop")
non_ent <- sample_ent |> transmute(sampleid, taxlev = "not Enterococcus",
                                   relab = 1 - ent_total)

keep <- sample_ent |> filter(ent_total >= 0.005) |> pull(sampleid)   # >= 0.5%
bar_df <- bind_rows(ent_long, non_ent) |> filter(sampleid %in% keep)

# sample ordering. "olo" replicates the published panel exactly: the original
# computed the order from a SPECIES-level barplot over ALL Enterococcus-containing
# samples and reused it for the ASV barplot. So we seriate the 16S Enterococcus
# *species* composition (878 samples) with optimal leaf ordering on ward.D2 /
# Bray-Curtis -- microViz comp_barplot's default -- then keep that order on the
# >= 0.5% subset. Other methods seriate the ASV composition of the kept samples.
comp_mat <- bar_df |>
  pivot_wider(names_from = taxlev, values_from = relab, values_fill = 0) |>
  column_to_rownames("sampleid")
order_samples <- function(method) {
  if (method == "olo") {
    # species-level composition = aggregate the ASV table to 16S species
    spmat <- e16 |>
      group_by(sampleid, species) |> summarise(relab = sum(relab), .groups = "drop") |>
      pivot_wider(names_from = species, values_from = relab, values_fill = 0) |>
      column_to_rownames("sampleid")
    ord <- seriation::seriate(vegdist(spmat, "bray"), method = "OLO_ward")
    return(rownames(spmat)[seriation::get_order(ord)])
  }
  d <- vegdist(comp_mat, "bray")
  switch(method,
    pcoa1   = rownames(comp_mat)[order(cmdscale(d, k = 1)[, 1])],
    hclust  = rownames(comp_mat)[hclust(d)$order],
    entdesc = sample_ent |> filter(sampleid %in% keep) |>
                arrange(desc(ent_total)) |> pull(sampleid),
    asv1    = comp_mat |> rownames_to_column("s") |> arrange(desc(asv_1)) |> pull(s))
}
sample_order <- order_samples(ORDER)
sample_order <- sample_order[sample_order %in% keep]   # restrict to plotted samples
if (tolower(Sys.getenv("BAR_REV", "false")) %in% c("true", "1")) sample_order <- rev(sample_order)

tax_levels <- c(top10, "Other Enterococcus", "not Enterococcus")
bar_df <- bar_df |>
  mutate(sampleid = factor(sampleid, levels = sample_order),
         taxlev = factor(taxlev, levels = tax_levels))

panel_c <- ggplot(bar_df, aes(sampleid, relab, fill = taxlev)) +
  geom_col(width = 1, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = asv_pal, drop = FALSE, name = "16S Enterococcus ASVs") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(y = "Abundance",
       x = sprintf("%d fecal samples with Enterococcus relative abundance >= 0.5%% by 16S sequencing",
                   length(keep)),
       title = "Stacked barplot of taxonomic composition by 16S sequencing") +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2,
                             byrow = TRUE, keywidth = unit(4, "mm"),
                             keyheight = unit(4, "mm"))) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 11),
        plot.title.position = "plot")

ggsave(out("E7c_asv_distribution.pdf"), panel_c, width = 9, height = 4.5)
message("wrote results/E7c_asv_distribution.pdf  (", length(keep),
        " samples, order = ", ORDER, ")")

# =====================================================================
# E7d: per-ASV NCBI RefSeq species matches, MGX prevalence + abundance
# =====================================================================
sp_abund <- read_csv(released("mgx_enterococcus_species_relab.csv"), show_col_types = FALSE)

# ASV -> amplicon-identical species (fixed lookup; see header), one assignment label
clusters <- tribble(
  ~asv,     ~assignment,         ~species,
  "asv_1",  "asv_1 faecium",     "Enterococcus_faecium",
  "asv_1",  "asv_1 faecium",     "Enterococcus_durans",
  "asv_1",  "asv_1 faecium",     "Enterococcus_hirae",
  "asv_1",  "asv_1 faecium",     "Enterococcus_mundtii",
  "asv_33", "asv_33 faecalis",   "Enterococcus_faecalis",
  "asv_33", "asv_33 faecalis",   "Enterococcus_dispar",
  "asv_67", "asv_67 gallinarum", "Enterococcus_gallinarum",
  "asv_67", "asv_67 gallinarum", "Enterococcus_casseliflavus",
  "asv_12", "asv_12 avium",      "Enterococcus_avium",
  "asv_12", "asv_12 avium",      "Enterococcus_gilvus",
  "asv_12", "asv_12 avium",      "Enterococcus_raffinosus")

prevalence <- sp_abund |> group_by(species) |>
  summarise(prev = sum(relab > 0), .groups = "drop")

row_meta <- clusters |>
  left_join(prevalence, by = "species") |>
  mutate(prev = replace_na(prev, 0),
         sp_short = str_remove(species, "Enterococcus_"),
         asv = factor(asv, levels = c("asv_1", "asv_33", "asv_67", "asv_12")),
         assignment = factor(assignment, levels = unique(assignment))) |>
  arrange(asv, desc(prev)) |>
  mutate(row_label = sprintf("italic('%s')~'  %d/%d'", sp_short, prev, N_MGX))
row_levels <- rev(row_meta$row_label)

points_df <- sp_abund |> filter(relab > 0) |>
  inner_join(row_meta, by = "species") |>
  mutate(row_label = factor(row_label, levels = row_levels))
meds <- points_df |> group_by(assignment, row_label) |>
  summarise(med = median(relab), .groups = "drop")

panel_d <- ggplot(points_df, aes(relab, row_label)) +
  geom_boxplot(outlier.colour = NA, fill = "gray85", colour = NA, width = 0.6) +
  geom_jitter(height = 0.2, width = 0, size = 1.1, alpha = 0.25, colour = "gray20") +
  geom_point(data = meds, aes(x = med, y = row_label), colour = "red",
             shape = 124, size = 5) +
  facet_grid(assignment ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1.0), expand = c(0.01, 0)) +
  labs(x = "MGX species abundance", y = NULL,
       caption = "grey box = IQR; red line = median; points = samples with nonzero abundance") +
  theme_classic(base_size = 11) +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "gray95", colour = NA),
        strip.text.y.left = element_text(angle = 0, face = "italic", size = 9),
        panel.spacing = unit(2, "mm"),
        plot.title = element_text(hjust = 0.5, size = 11),
        plot.caption = element_text(hjust = 0.5, size = 8),
        axis.text.y = element_text(size = 9))

ggsave(out("E7d_species_table.pdf"), panel_d, width = 9, height = 6)
message("wrote results/E7d_species_table.pdf")
