# Figure 1 n,o: diet vs microbiome beta diversity.
#
# Refactor of R06_beta_diversity_diet_and_microbiome. For each stool sample we take
# the prior-2-day average diet (food-code dehydrated weight), and compare diet
# distance to microbiome distance within each patient, relative to that patient's
# earliest sample.
#
#   F1n  two PCoA scatters side by side -- Diet PCoA (food-code unweighted UniFrac)
#        and Microbiome PCoA (genus Bray-Curtis) -- with one illustrative patient
#        (P52) highlighted (earliest timepoint = triangle, later = circles).
#   F1o  within-patient microbiome distance-from-earliest vs diet distance-from-
#        earliest, with the linear mixed model fit (stool_dist ~ diet_dist + (1|pid)).
#
# Diet UniFrac runs in one QIIME2 Docker container (same pattern as the procrustes
# script); the microbiome Bray-Curtis PCoA is done in R (vegan), as in R06.
#
# Substitution: R06 read 022_ALL173_stool_samples_genus_counts.csv, which is not
# shipped; the genus relative abundance is rebuilt from the released per-ASV table
# (171_quality_asv_relab_pident97_genus.csv), as in F4a/E7b. So the microbiome PCoA
# variance and the F1o fit may differ slightly from the published panel. FLAGGED.
#
# RUN_QIIME=false reuses the cached diet ordination/distance; QIIME2_IMAGE picks tag.

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggpubr)
  library(lmerTest)
  library(patchwork)
})

ROOT      <- here::here()
released  <- function(f) file.path(ROOT, "released_data", f)
WORK      <- file.path(ROOT, "intermediate_data", "R06_diet_beta")
TREE_NWK  <- released("output_food_tree_datatree.newick")
QIIME2_IMAGE <- Sys.getenv("QIIME2_IMAGE", unset = "quay.io/qiime2/qiime2:2026.4")
RUN_QIIME    <- tolower(Sys.getenv("RUN_QIIME", unset = "true")) %in% c("true", "1", "yes")
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ROOT, "results"), showWarnings = FALSE)

dtb  <- read_csv(released("152_combined_DTB.csv"),  show_col_types = FALSE)
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)

# ---- 1. prior-2-day average food-code matrix (rows = Food_code, cols = sample) -
food_daily <- dtb |> group_by(pid, fdrt, Food_code) |>
  summarise(w = sum(dehydrated_weight), .groups = "drop")
p2d_fc <- meta |> select(pid, sdrt, sampleid) |>
  mutate(ws = sdrt - 2, we = sdrt - 1) |>
  left_join(food_daily, by = join_by(pid, ws <= fdrt, we >= fdrt)) |>
  filter(!is.na(Food_code)) |>
  group_by(sampleid, Food_code) |> summarise(ave_fc = sum(w) / 2, .groups = "drop")

fc_mat <- p2d_fc |>
  pivot_wider(names_from = sampleid, values_from = ave_fc, values_fill = 0) |>
  rename(`#OTU ID` = Food_code)
fc_mat <- fc_mat[rowSums(as.matrix(fc_mat[, -1])) > 0, ]
write_tsv(fc_mat, file.path(WORK, "R06_p2d_fc_df.tsv"))
message(sprintf("diet matrix: %d food codes x %d samples", nrow(fc_mat), ncol(fc_mat) - 1))

# ---- 2. QIIME2: diet unweighted UniFrac -> distance matrix + PCoA --------------
run_qiime <- function() {
  rel <- function(p) sub(paste0("^", ROOT, "/?"), "", p)
  base_c <- file.path("/data", rel(WORK)); tree_c <- file.path("/data", rel(TREE_NWK))
  script <- sprintf('
set -euo pipefail
pip uninstall -y q2-composition >/dev/null 2>&1 || true
BASE=%s
biom convert -i "$BASE/R06_p2d_fc_df.tsv" -o "$BASE/fc.biom" --to-hdf5 --table-type="Table"
qiime tools import --input-path "$BASE/fc.biom" --output-path "$BASE/fc.qza" --type "FeatureTable[Frequency]"
[ -f "$BASE/food_tree.qza" ] || qiime tools import --input-path %s --output-path "$BASE/food_tree.qza" --type "Phylogeny[Rooted]"
qiime diversity beta-phylogenetic --i-table "$BASE/fc.qza" --i-phylogeny "$BASE/food_tree.qza" --p-metric unweighted_unifrac --o-distance-matrix "$BASE/diet_dm.qza"
qiime diversity pcoa --i-distance-matrix "$BASE/diet_dm.qza" --o-pcoa "$BASE/diet_pcoa.qza"
rm -rf "$BASE/diet_dm_export" "$BASE/diet_pcoa_export"
qiime tools export --input-path "$BASE/diet_dm.qza"   --output-path "$BASE/diet_dm_export"
qiime tools export --input-path "$BASE/diet_pcoa.qza" --output-path "$BASE/diet_pcoa_export"
', base_c, tree_c)
  cmd <- sprintf('docker run --rm --platform linux/amd64 -v %s:/data -w /data %s bash -c %s',
                 shQuote(ROOT), shQuote(QIIME2_IMAGE), shQuote(script))
  message("running QIIME2 diet UniFrac + PCoA ...")
  if (system(cmd) != 0) stop("QIIME2 container exited non-zero")
}
dm_tsv  <- file.path(WORK, "diet_dm_export", "distance-matrix.tsv")
ord_txt <- file.path(WORK, "diet_pcoa_export", "ordination.txt")
if (RUN_QIIME || !file.exists(ord_txt)) run_qiime() else message("reusing cached diet ordination")

# ---- helpers to parse the exported ordination ---------------------------------
read_site <- function(f) {
  L <- readLines(f); h <- grep("^Site\t", L)[1]
  n <- as.integer(strsplit(L[h], "\t")[[1]][2]); b <- strsplit(L[(h + 1):(h + n)], "\t")
  tibble(sampleid = vapply(b, `[`, character(1), 1),
         PC1 = as.numeric(vapply(b, `[`, character(1), 2)),
         PC2 = as.numeric(vapply(b, `[`, character(1), 3)))
}
read_prop <- function(f) {
  L <- readLines(f); i <- grep("^Proportion explained", L)[1]
  as.numeric(strsplit(L[i + 1], "\t")[[1]])
}

# within-patient distance from the earliest sample (R06 logic)
dist_from_earliest <- function(long) {
  long <- long |>
    left_join(meta |> select(pid1 = pid, sampleid1 = sampleid, sdrt1 = sdrt), by = "sampleid1") |>
    left_join(meta |> select(pid2 = pid, sampleid2 = sampleid, sdrt2 = sdrt), by = "sampleid2") |>
    filter(pid1 == pid2, sampleid1 != sampleid2) |> rename(pid = pid1) |> select(-pid2)
  min_sdrt <- long |> group_by(pid) |> summarise(min_sdrt = min(c(sdrt1, sdrt2)), .groups = "drop")
  long |> left_join(min_sdrt, by = "pid") |>
    filter(sdrt1 == min_sdrt, sdrt2 != min_sdrt)
}

# ---- 3. diet: PCoA + distances -----------------------------------------------
diet_var <- round(read_prop(ord_txt) * 100, 1)
diet_pcoa <- read_site(ord_txt) |> inner_join(meta |> select(sampleid, pid, sdrt), by = "sampleid")

diet_long <- read_tsv(dm_tsv, show_col_types = FALSE) |> rename(sampleid1 = 1) |>
  pivot_longer(-sampleid1, names_to = "sampleid2", values_to = "distance")
diet_dist <- dist_from_earliest(diet_long) |> rename(diet_dist = distance)

# ---- 4. microbiome: genus Bray-Curtis PCoA + distances (R/vegan) --------------
genus_relab <- read_csv(released("171_quality_asv_relab_pident97_genus.csv"), show_col_types = FALSE) |>
  filter(!is.na(genus), sampleid %in% meta$sampleid) |>
  group_by(sampleid, genus) |> summarise(relab = sum(count_relative), .groups = "drop") |>
  pivot_wider(names_from = genus, values_from = relab, values_fill = 0) |>
  column_to_rownames("sampleid")

bray <- vegdist(genus_relab, "bray")
eig <- cmdscale(bray, eig = TRUE)$eig
stool_var <- signif(eig / sum(eig), 3) * 100
stool_pcoa <- cmdscale(bray, k = 2) |> as.data.frame() |> rownames_to_column("sampleid") |>
  rename(PC1 = V1, PC2 = V2) |> inner_join(meta |> select(sampleid, pid, sdrt), by = "sampleid")

stool_long <- as.matrix(bray) |> as.data.frame() |> rownames_to_column("sampleid1") |>
  pivot_longer(-sampleid1, names_to = "sampleid2", values_to = "distance")
stool_dist <- dist_from_earliest(stool_long) |> rename(stool_dist = distance)

# ---- 5. F1o: microbiome vs diet distance, linear mixed model ------------------
joined <- stool_dist |>
  inner_join(diet_dist |> select(pid, sampleid1, sampleid2, diet_dist),
             by = c("pid", "sampleid1", "sampleid2"))
model <- lmer(stool_dist ~ diet_dist + (1 | pid), data = joined)
co <- summary(model)$coefficients
b0 <- co[1, 1]; b1 <- co[2, 1]; pval <- co[2, 5]
message(sprintf("F1o lmer: stool = %.3f + %.3f*diet,  p = %.2g  (n=%d pairs, %d patients)",
                b0, b1, pval, nrow(joined), n_distinct(joined$pid)))

# Plot the PARTIAL RESIDUALS (microbiome distance with the per-patient random
# intercept removed): residuals + fixed-effect prediction. This is what R06 plots;
# it disperses the cloud around the fixed-effect line and is not capped at [0,1],
# unlike the raw Bray-Curtis distance (which piles up near 1).
joined <- joined |>
  mutate(partial = residuals(model) + predict(model, re.form = NA))

f1o <- ggplot(joined, aes(diet_dist, partial)) +
  geom_point(colour = "#7B68A6", alpha = 0.35, size = 1.4) +
  geom_abline(intercept = b0, slope = b1, linewidth = 1) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4, size = 3.4,
           label = sprintf("y = %.3f + %.3fx\np %s", b0, b1,
                           ifelse(pval < 0.001, "< 0.001", sprintf("= %.3f", pval)))) +
  labs(x = "Diet Distance", y = expression("microbiome " * beta * "-diversity")) +
  theme_classic() + theme(aspect.ratio = 1)

# ---- 6. F1n: the two PCoA scatters with patient P52 highlighted ----------------
selected <- "P52"
earliest_sdrt <- meta |> filter(pid == selected) |> slice_min(sdrt, n = 1) |> pull(sdrt)
tag <- function(df) df |>
  mutate(grp = case_when(pid == selected & sdrt == earliest_sdrt ~ "earliest",
                         pid == selected ~ "later", TRUE ~ "other")) |>
  arrange(grp != "other")

pcoa_panel <- function(df, bgcol, title, xlab, ylab) {
  df <- tag(df)
  ggplot(mapping = aes(PC1, PC2)) +
    geom_point(data = filter(df, grp == "other"), colour = bgcol, alpha = 0.45, size = 1.1) +
    geom_point(data = filter(df, grp == "later"),    colour = "red", shape = 16, size = 2.6) +
    geom_point(data = filter(df, grp == "earliest"), colour = "red", shape = 17, size = 3.4) +
    labs(title = title, x = xlab, y = ylab) +
    theme_classic() +
    theme(aspect.ratio = 1, plot.title = element_text(size = 10, hjust = 0.5))
}
p_diet  <- pcoa_panel(diet_pcoa, "gray35", "Diet PCoA",
                      sprintf("PCo-1 [%s%%]", diet_var[1]), sprintf("PCo-2 [%s%%]", diet_var[2]))
p_stool <- pcoa_panel(stool_pcoa, "#4F6FB0", "Microbiome PCoA",
                      sprintf("PCo-1 [%s%%]", stool_var[1]), sprintf("PCo-2 [%s%%]", stool_var[2]))

# shared legend strip at the bottom
leg_df <- tibble(x = 1:4, lab = c("diet-data days", "fecal samples",
                                  "P52 earliest", "P52 later"),
                 col = c("gray35", "#4F6FB0", "red", "red"), sh = c(16, 16, 17, 16))
legend <- ggplot(leg_df, aes(x, 1)) +
  geom_point(aes(colour = I(col), shape = I(sh)), size = 3) +
  geom_text(aes(label = lab), vjust = 2.4, size = 3) +
  scale_x_continuous(limits = c(0.5, 4.6)) + ylim(0.6, 1.1) + theme_void()

f1n <- (p_diet | p_stool)
ggsave(file.path(ROOT, "results", "F1n_diet_microbiome_pcoa.pdf"),
       f1n / legend + plot_layout(heights = c(1, 0.18)), width = 200, height = 120, units = "mm")
ggsave(file.path(ROOT, "results", "F1o_diet_stool_distance.pdf"), f1o,
       width = 95, height = 95, units = "mm")
message("wrote results/F1n_diet_microbiome_pcoa.pdf and results/F1o_diet_stool_distance.pdf")
message(sprintf("\nDiet PCoA var: PCo1 %.1f%%  PCo2 %.1f%%   |   Microbiome PCoA var: PCo1 %.1f%%  PCo2 %.1f%%",
                diet_var[1], diet_var[2], stool_var[1], stool_var[2]))
