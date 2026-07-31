# Figure 2a: how many prior days of diet to correlate with each stool sample.
#
# For each window W (1..5 prior days), diet and stool are ordinated and matched by a
# symmetric Procrustes fit:
#   food   food-code dehydrated weight, daily-averaged over the W prior days
#   macro  the five macronutrients, daily-averaged over the W prior days
#   stool  stool composition collapsed to genus (window-independent)
# QIIME2 turns each into a PCoA ordination (food -> unweighted UniFrac on the food
# tree; macro/stool -> Bray-Curtis).
#
# The panel plots, on a shared x (diet exposure days):
#   left axis (black)   Delta Procrustes correlation, the MARGINAL GAIN in the
#                       symmetric Procrustes correlation sqrt(1 - M^2) between
#                       successive windows, x 10^-2
#   right axis (grey)   the Procrustes correlation itself
# for FNDDS food groups (solid, circle) and macronutrients (dashed, square). A
# delete-one-patient jackknife tests whether the gain from one one-day extension
# differs from the next; the 1->2 food-group step is the largest and significant.
#
# Cohort: the 801 stool samples with all five prior days recorded (a food entry, or
# a documented 072 zero-eating day), minus those whose most recent prior day is a
# zero-eating day (empty food composition, so no UniFrac distance), leaving a fixed
# set of 785 samples / 143 patients scored identically across all five windows.
#
# The QIIME work runs once inside a single container and caches its exported
# ordinations under intermediate_data/006_paired_for_procrustes<suffix>/ (idempotent;
# RUN_QIIME=false reuses them). The jackknife caches to
# intermediate_data/006_jackknife_pvalues<suffix>.csv (RUN_JACKKNIFE=false reuses it).
# Panel: results/F2a_procrustes_delta<suffix>.pdf.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages(library(vegan))

# INCLUDE_ZERO_DAYS (default true): treat documented zero-eating days (072) as valid
# zero-intake data rather than missing, matching the E1a cohort logic and the
# manuscript Methods. A prior day counts as "covered" if it has a dietary entry OR a
# 072 zero-eating record; this is the 801-then-785 cohort above. Set false to build
# the entry-only cohort instead (751), which writes to unsuffixed cache / file names.
INCLUDE_ZERO_DAYS <- tolower(Sys.getenv("INCLUDE_ZERO_DAYS", "true")) %in% c("true","1","yes")
suffix <- if (INCLUDE_ZERO_DAYS) "_zerodays" else ""

PAIR_ROOT <- cache_path(paste0("006_paired_for_procrustes", suffix))
TREE_NWK  <- released("output_food_tree_datatree.newick")   # Food_code leaves, rooted
QIIME2_IMAGE <- Sys.getenv("QIIME2_IMAGE", unset = "quay.io/qiime2/qiime2:2026.4")
RUN_QIIME    <- tolower(Sys.getenv("RUN_QIIME", unset = "true")) %in% c("true", "1", "yes")

dir.create(PAIR_ROOT, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(TREE_NWK))

MACROS <- c("Protein_g", "Fat_g", "Carbohydrates_g", "Fibers_g", "Sugars_g")
Ws <- 1:5

# ---- 1. inputs -------------------------------------------------------------
dtb  <- read_csv(released("152_combined_DTB.csv"),  show_col_types = FALSE)
meta <- read_csv(released("153_combined_META.csv"), show_col_types = FALSE)
asv  <- read_csv(released("45_quality_asv_relab_pident97_genus.csv"),
                 show_col_types = FALSE)

# ---- 2. samples with a food record on all 5 prior days ---------------------
# A prior day is "covered" if it has a dietary entry, or (when INCLUDE_ZERO_DAYS)
# is a documented zero-eating day from 072 (real zero intake, not missing data).
diet_days <- dtb |> distinct(pid, fdrt)
if (INCLUDE_ZERO_DAYS) {
  zero_days <- read_csv(released("072_total_patients_zero_eating_days_pid.csv"),
                        show_col_types = FALSE) |> distinct(pid, fdrt)
  diet_days <- bind_rows(diet_days, zero_days) |> distinct(pid, fdrt)
}
diet_days <- diet_days |> mutate(has = TRUE)

prior_ok <- meta |>
  select(pid, sdrt, sampleid) |>
  tidyr::crossing(k = 1:5) |>
  mutate(fdrt = sdrt - k) |>
  left_join(diet_days, by = c("pid", "fdrt")) |>
  group_by(sampleid) |>
  summarise(n_ok = sum(!is.na(has)), .groups = "drop") |>
  filter(n_ok == 5)

meta_keep    <- meta |> semi_join(prior_ok, by = "sampleid")
sids_ordered <- sort(unique(meta_keep$sampleid))
message("Samples kept: ", nrow(meta_keep), " / ", nrow(meta))

# ---- 3. build feature tables and write the QIIME TSVs ----------------------
long_diet <- meta_keep |>
  select(pid, sdrt, sampleid) |>
  tidyr::crossing(offset = 1:5) |>
  mutate(fdrt = sdrt - offset) |>
  inner_join(
    dtb |> group_by(pid, fdrt, Food_code) |>
      summarise(dehydrated_weight = sum(dehydrated_weight), .groups = "drop"),
    by = c("pid", "fdrt"), relationship = "many-to-many"
  )

long_macro <- meta_keep |>
  select(pid, sdrt, sampleid) |>
  tidyr::crossing(offset = 1:5) |>
  mutate(fdrt = sdrt - offset) |>
  inner_join(
    dtb |> group_by(pid, fdrt) |>
      summarise(across(all_of(MACROS), sum), .groups = "drop"),
    by = c("pid", "fdrt")
  )

# Pad to every kept sample column in fixed order, drop all-zero feature rows AND
# all-zero sample columns (a sample with no diet in this window has an undefined
# UniFrac/Bray distance, so it cannot be ordinated; the Procrustes intersects
# common samples). Rename the id column to the '#OTU ID' header biom convert wants.
finalize_table <- function(tbl, id_col) {
  for (s in setdiff(sids_ordered, names(tbl)[-1])) tbl[[s]] <- 0
  tbl <- tbl[, c(id_col, sids_ordered)]
  M <- as.matrix(tbl[, -1])
  tbl <- tbl[rowSums(M) > 0, c(TRUE, colSums(M) > 0), drop = FALSE]
  names(tbl)[1] <- "#OTU ID"
  tbl
}

# Stool genus relative abundance (window-independent, built once).
genus_tbl <- asv |>
  filter(!is.na(genus), sampleid %in% sids_ordered) |>
  group_by(sampleid, genus) |>
  summarise(count_relative = sum(count_relative), .groups = "drop") |>
  pivot_wider(names_from = sampleid, values_from = count_relative, values_fill = 0) |>
  finalize_table("genus")

for (W in Ws) {
  P <- file.path(PAIR_ROOT, paste0("p", W, "day")); T <- paste0("d", W)
  dir.create(P, recursive = TRUE, showWarnings = FALSE)

  food <- long_diet |>
    filter(offset <= W) |>
    group_by(sampleid, Food_code) |>
    summarise(daily_avg = sum(dehydrated_weight) / W, .groups = "drop") |>
    pivot_wider(names_from = sampleid, values_from = daily_avg, values_fill = 0) |>
    finalize_table("Food_code")

  macro <- long_macro |>
    filter(offset <= W) |>
    group_by(sampleid) |>
    summarise(across(all_of(MACROS), ~ sum(.x) / W), .groups = "drop") |>
    pivot_longer(all_of(MACROS), names_to = "macro", values_to = "daily_avg") |>
    pivot_wider(names_from = sampleid, values_from = daily_avg, values_fill = 0) |>
    finalize_table("macro")

  write_tsv(food,      file.path(P, paste0(T, "_food.tsv")))
  write_tsv(macro,     file.path(P, paste0(T, "_macro.tsv")))
  write_tsv(genus_tbl, file.path(P, paste0(T, "_stool_genus.tsv")))
  message(sprintf("W=%d  food=%d x %d  macro=%d x %d  stool_genus=%d",
                  W, nrow(food), ncol(food) - 1, nrow(macro), ncol(macro) - 1,
                  nrow(genus_tbl)))
}

# ---- 4. QIIME2: one container, whole W loop --------------------------------
# tsv -> biom -> qza -> beta -> PCoA -> export. The food tree is imported once.
run_qiime <- function() {
  root <- here::here()
  rel  <- function(p) sub(paste0("^", root, "/?"), "", p)        # host path -> /data path
  base_c <- file.path("/data", rel(PAIR_ROOT))
  tree_c <- file.path("/data", rel(TREE_NWK))

  script <- sprintf('
set -euo pipefail
# The 2026.4 image ships a q2-composition that aborts the plugin manager trying to
# import the R package phyloseq (not installed). We do not use composition; drop it
# so qiime commands load. The container is --rm, so this is ephemeral.
pip uninstall -y q2-composition >/dev/null 2>&1 || true
BASE=%s
[ -f "$BASE/food_tree.qza" ] || qiime tools import --input-path %s --output-path "$BASE/food_tree.qza" --type "Phylogeny[Rooted]"

# tsv -> biom -> qza -> beta -> PCoA -> export, skipping anything already exported.
# $1 = table stem (no ext), $2 = export dir, $3 = metric, $4 = tree qza (phylo only)
ensure_pcoa () {
  if [ -f "$2/ordination.txt" ]; then echo "skip $2"; return; fi
  biom convert -i "$1.tsv" -o "$1.biom" --to-hdf5 --table-type="Table"
  qiime tools import --input-path "$1.biom" --output-path "$1.qza" --type "FeatureTable[Frequency]"
  if [ -n "${4:-}" ]; then
    qiime diversity beta-phylogenetic --i-table "$1.qza" --i-phylogeny "$4" --p-metric "$3" --o-distance-matrix "$1_dm.qza"
  else
    qiime diversity beta --i-table "$1.qza" --p-metric "$3" --o-distance-matrix "$1_dm.qza"
  fi
  qiime diversity pcoa --i-distance-matrix "$1_dm.qza" --o-pcoa "$1_dm_pcoa.qza"
  qiime tools export --input-path "$1_dm_pcoa.qza" --output-path "$2"
}

for W in 1 2 3 4 5; do
  P="$BASE/p${W}day"; T="d${W}"
  echo "===== W=${W} ====="
  ensure_pcoa "$P/${T}_food"        "$P/${T}_food_pcoa"        unweighted_unifrac "$BASE/food_tree.qza"
  ensure_pcoa "$P/${T}_macro"       "$P/${T}_macro_pcoa"       braycurtis
  ensure_pcoa "$P/${T}_stool_genus" "$P/${T}_stool_genus_pcoa" braycurtis
done
', base_c, tree_c)

  cmd <- sprintf(
    'docker run --rm --platform linux/amd64 -v %s:/data -w /data %s bash -c %s',
    shQuote(root), shQuote(QIIME2_IMAGE), shQuote(script))
  message("running QIIME2 in one container ...")
  status <- system(cmd)
  if (status != 0) stop("QIIME2 container exited with status ", status)
}

ordo_dir <- function(W, kind) file.path(PAIR_ROOT, paste0("p", W, "day"),
                                        sprintf("d%d_%s_pcoa", W, kind))
KINDS <- c("food", "macro", "stool_genus")
have_all <- all(file.exists(file.path(
  as.vector(outer(Ws, KINDS, function(w, k) mapply(ordo_dir, w, k))), "ordination.txt")))

if (RUN_QIIME || !have_all) run_qiime() else message("reusing cached ordinations")

# ---- 5. read ordinations; fixed cohort; symmetric Procrustes correlation ---
# Parse the skbio ordination.txt 'Site' block into a sample x PC matrix.
read_pcoa_matrix <- function(dir) {
  lines <- readLines(file.path(dir, "ordination.txt"))
  h <- grep("^Site\t", lines)[1]
  n <- as.integer(strsplit(lines[h], "\t")[[1]][2])
  block <- strsplit(lines[(h + 1):(h + n)], "\t")
  mat <- do.call(rbind, lapply(block, function(r) as.numeric(r[-1])))
  rownames(mat) <- vapply(block, `[`, character(1), 1)
  mat
}

# Load all 15 ordinations (5 windows x food/macro/stool) once.
ord <- list()
for (kind in KINDS) for (W in Ws) ord[[paste(kind, W)]] <- read_pcoa_matrix(ordo_dir(W, kind))

# Fixed analysis cohort: samples present in EVERY ordination, so all five windows
# are scored on the same set. The day-1 food ordination already drops the samples
# whose most recent prior day has no food composition, so the intersection is 785.
ids <- Reduce(intersect, lapply(ord, rownames))
samp_pid <- meta |> filter(sampleid %in% ids) |> distinct(sampleid, pid)
message(sprintf("Analysis cohort: %d samples, %d patients",
                length(ids), dplyr::n_distinct(samp_pid$pid)))

diet_kinds <- c("Food group based" = "food", "Macronutrient based" = "macro")

# Symmetric Procrustes correlation sqrt(1 - M^2) per window over a set of samples.
# symmetric = TRUE scales each ordination to unit sum of squares, so M^2 in [0,1].
corr_windows <- function(kind, keep) {
  vapply(Ws, function(W) {
    d <- ord[[paste(kind, W)]]; s <- ord[[paste("stool_genus", W)]]
    cc <- intersect(keep, intersect(rownames(d), rownames(s)))
    sqrt(1 - vegan::procrustes(d[cc, ], s[cc, ], symmetric = TRUE)$ss)
  }, numeric(1))
}

# Correlation per window and the marginal gain (Delta) between successive windows.
corr_tbl <- purrr::imap_dfr(diet_kinds, function(kind, pr)
  tibble(pair = pr, W = Ws, corr = corr_windows(kind, ids))) |>
  group_by(pair) |> arrange(W, .by_group = TRUE) |>
  mutate(gain = corr - lag(corr)) |> ungroup()

write_csv(corr_tbl, cache_path(paste0("006_procrustes_delta", suffix, ".csv")))
message("\nProcrustes correlation and marginal gain (x10^-2) by window:")
print(corr_tbl |> mutate(gain_x100 = round(gain * 100, 3)), n = Inf)

# ---- 6. delete-one-patient jackknife for the marginal-gain comparisons -----
# Statistic per diet: gain(W -> W+1) - gain(W+1 -> W+2), i.e. whether one one-day
# extension adds more than the next. Two-sided z-test with a delete-one-patient
# jackknife standard error. Cached; RUN_JACKKNIFE=false reuses the cache.
gains_of    <- function(cv) cv[-1] - cv[-length(cv)]        # 1->2, 2->3, 3->4, 4->5
consec_diff <- function(g)  g[-length(g)] - g[-1]           # (1->2)-(2->3), ...
JACK_CACHE  <- cache_path(paste0("006_jackknife_pvalues", suffix, ".csv"))
RUN_JACK    <- tolower(Sys.getenv("RUN_JACKKNIFE", "true")) %in% c("true", "1", "yes")

jackknife_p <- function() {
  pats  <- unique(samp_pid$pid); n <- length(pats)
  steps <- c("(1->2)-(2->3)", "(2->3)-(3->4)", "(3->4)-(4->5)")
  purrr::imap_dfr(diet_kinds, function(kind, pr) {
    full <- consec_diff(gains_of(corr_windows(kind, ids)))
    J <- vapply(pats, function(p)
      consec_diff(gains_of(corr_windows(kind, samp_pid$sampleid[samp_pid$pid != p]))),
      numeric(length(steps)))
    Jt <- t(J)                                              # n patients x 3 steps
    se <- sqrt((n - 1) / n * colSums(sweep(Jt, 2, colMeans(Jt))^2))
    z  <- full / se
    tibble(pair = pr, step = steps, d = full, se = se, z = z, p = 2 * pnorm(-abs(z)))
  })
}

jack <- if (!RUN_JACK && file.exists(JACK_CACHE)) {
  message("reusing cached jackknife p-values")
  read_csv(JACK_CACHE, show_col_types = FALSE)
} else {
  message(sprintf("running delete-one-patient jackknife (%d patients x 5 windows x 2 diets) ...",
                  dplyr::n_distinct(samp_pid$pid)))
  j <- jackknife_p(); write_csv(j, JACK_CACHE); j
}
message("\nMarginal-gain comparisons (jackknife):")
print(jack, n = Inf)

# p value annotated on the panel: food-group 1->2 vs 2->3 gain comparison.
p_food_12 <- jack |> filter(pair == "Food group based", step == "(1->2)-(2->3)") |> pull(p)

# ---- 7. F2a: marginal gain (left, black) + correlation (right, grey) -------
# Left axis is the marginal gain x10^-2; the right axis (correlation, 0.2..1.0) is
# mapped onto the left coordinate by L = K * corr, so K = 2.5 puts corr = 1.0 at the
# top of the 0..2.5 gain range. Delta points sit at the window midpoint (1.5..4.5),
# dodged slightly by diet; correlation points at the integer days.
K <- 2.5
delta <- corr_tbl |> filter(!is.na(gain)) |>
  mutate(x = W - 0.5 + ifelse(pair == "Food group based", -0.06, 0.06), y = gain * 100)
grey  <- corr_tbl |> mutate(y = K * corr)

lt <- c("Food group based" = "solid", "Macronutrient based" = "dashed")
sh <- c("Food group based" = 16, "Macronutrient based" = 15)

# Conventional threshold annotation (manuscript style, capital italic P) rather than
# the exact jackknife value: the tightest 10^-n / 0.05 threshold the p value clears.
# the thresholds are quoted inside the expression so plotmath prints them literally
# (0.0001) rather than reformatting the number to scientific notation (1e-04).
plab <- if (p_food_12 < 1e-4) 'italic(P) < "0.0001"' else
        if (p_food_12 < 1e-3) 'italic(P) < "0.001"'  else
        if (p_food_12 < 1e-2) 'italic(P) < "0.01"'   else
        if (p_food_12 < 0.05) 'italic(P) < "0.05"'   else
        sprintf('italic(P) == "%.2f"', p_food_12)

f2a <- ggplot() +
  # grey: the correlation itself, on the right axis
  geom_line(data = grey, aes(W, y, group = pair, linetype = pair),
            colour = "grey65", linewidth = 0.5) +
  geom_point(data = grey, aes(W, y, shape = pair), colour = "grey65", size = 1.6) +
  # black: the marginal gain, on the left axis
  geom_line(data = delta, aes(x, y, group = pair, linetype = pair), linewidth = 0.9) +
  geom_point(data = delta, aes(x, y, shape = pair), size = 2.6) +
  # bracket + p value over the food-group 1->2 vs 2->3 comparison
  annotate("segment", x = 1.44, xend = 2.44, y = 2.85, yend = 2.85, linewidth = 0.4) +
  annotate("segment", x = 1.44, xend = 1.44, y = 2.78, yend = 2.85, linewidth = 0.4) +
  annotate("segment", x = 2.44, xend = 2.44, y = 2.78, yend = 2.85, linewidth = 0.4) +
  annotate("text", x = 1.94, y = 2.98, label = plab, size = 4, parse = TRUE) +
  scale_shape_manual(values = sh) +
  scale_linetype_manual(values = lt) +
  scale_x_continuous(breaks = Ws, limits = c(0.8, 5.2)) +
  scale_y_continuous(
    name = expression(Delta ~ "Procrustes correlation (" * "×" * 10^-2 * ")"),
    breaks = 0:3, limits = c(0, 3.1),
    sec.axis = sec_axis(~ . / K, name = "Procrustes correlation",
                        breaks = c(0.2, 0.4, 0.6, 0.8, 1.0))) +
  labs(x = "diet exposure days", shape = NULL, linetype = NULL) +
  theme_classic() +
  theme(aspect.ratio = 1,
        plot.margin = margin(6, 6, 6, 12),
        legend.position = "inside", legend.position.inside = c(0.62, 0.82),
        legend.key.width = unit(1.1, "cm"),
        legend.background = element_blank(),
        axis.text  = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title_size),
        axis.title.y.right = element_text(colour = "grey55"),
        axis.text.y.right  = element_text(colour = "grey55"),
        axis.line.y.right  = element_line(colour = "grey55"),
        axis.ticks.y.right = element_line(colour = "grey55"))

# save_panel() uses the base pdf() device, which drops the plotmath Greek Delta in
# the y-axis title; cairo_pdf renders it (and the x10^-2 superscript) correctly.
f2a_path <- here::here("results", paste0("F2a_procrustes_delta", suffix, ".pdf"))
ggsave(f2a_path, f2a, width = 120, height = 120, units = "mm", device = cairo_pdf)
message("wrote ", f2a_path)
