# Figure 2a: how well diet composition predicts stool composition as the diet
# exposure window grows from 1 to 5 prior days, measured by Procrustes fit.
#
# Refactor of 006_procrustes.Rmd into the human family layout. For each window W
# (1..5 prior days), over the stool samples with a food record on all 5 prior days:
#   food   food-code dehydrated weight, daily-averaged over the W prior days
#   macro  the five macronutrients, daily-averaged over the W prior days
#   stool  stool composition collapsed to genus (window-independent)
# QIIME2 turns each into a PCoA ordination (food -> unweighted UniFrac on the food
# tree; macro/stool -> Bray-Curtis), and a Procrustes test scores how closely the
# diet ordination matches the stool ordination.
#
# Two Procrustes scorings are reported as SEPARATE panels (different file names):
#   asymmetric  vegan::procrustes(X = diet, Y = stool), default symmetric = FALSE;
#               SS on the raw ordination scale.
#   symmetric   vegan::procrustes(..., symmetric = TRUE); normalized M^2 in [0,1],
#               comparable across windows. The plotted metric is best-window-highest
#               (symmetric -> Procrustes correlation sqrt(1 - M^2); asymmetric ->
#               max(SS) - SS), so the tightest-fitting window is the top point.
#
# The QIIME work runs once inside a single container (one cold start, not ~75) and
# its exported ordinations cache under intermediate_data/006_paired_for_procrustes/.
# The container step is idempotent (skips any ordination already exported).
# Set RUN_QIIME=false to reuse cached ordinations; QIIME2_IMAGE picks the tag.

source(here::here("reproduce", "human", "_human_helpers.R"))
suppressPackageStartupMessages(library(vegan))

# INCLUDE_ZERO_DAYS: treat documented zero-eating days (072) as valid zero-intake
# data rather than missing. A prior day then counts as "covered" if it has either a
# dietary entry OR a 072 zero-eating record, matching the E1a cohort logic (where
# zero-eating days are not missing data). This enlarges the all-5-prior-days fixed
# cohort from 751 (entry-only) to 801, with zero-eating days contributing 0 intake
# to the windowed diet average. Outputs go to a separate cache / file names.
INCLUDE_ZERO_DAYS <- tolower(Sys.getenv("INCLUDE_ZERO_DAYS", "false")) %in% c("true","1","yes")
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

# ---- 5. read ordinations and run the Procrustes tests ----------------------
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

# Both Procrustes SS in one pass:
#   asymmetric = procrustes(X = diet, Y = stool)  (default symmetric = FALSE);
#                SS on the raw ordination scale.
#   symmetric  = procrustes(..., symmetric = TRUE); normalized M^2 in [0,1].
procrustes_test <- function(diet_dir, stool_dir) {
  d <- read_pcoa_matrix(diet_dir)
  s <- read_pcoa_matrix(stool_dir)
  common <- intersect(rownames(d), rownames(s))
  tibble(
    n = length(common),
    asymmetric = vegan::procrustes(X = d[common, ], Y = s[common, ])$ss,
    symmetric  = vegan::procrustes(X = d[common, ], Y = s[common, ], symmetric = TRUE)$ss
  )
}

diet_kinds <- c("Food group based" = "food", "Macronutrient based" = "macro")

procrustes_all <- tidyr::expand_grid(W = Ws, pair = names(diet_kinds)) |>
  mutate(pNd       = paste0("p", W, "d"),
         diet_dir  = ordo_dir(W, diet_kinds[pair]),
         stool_dir = ordo_dir(W, "stool_genus"),
         res       = map2(diet_dir, stool_dir, procrustes_test)) |>
  unnest(res) |>
  select(pNd, W, pair, n, asymmetric, symmetric) |>
  pivot_longer(c(asymmetric, symmetric), names_to = "method", values_to = "ss") |>
  group_by(pair, method) |>
  # best-window-highest metric: symmetric -> Procrustes correlation sqrt(1 - M^2)
  # (a proper fit in [0,1]); asymmetric SS is unbounded so use max(SS) - SS.
  mutate(score = if (first(method) == "symmetric") sqrt(1 - ss) else max(ss) - ss) |>
  ungroup()

write_csv(procrustes_all, cache_path(paste0("006_procrustes_scores", suffix, ".csv")))
message("\nProcrustes SS / score by method (genus stool):")
print(procrustes_all |> arrange(method, pair, W), n = Inf)

# Per-window available cohort: number of stool samples whose prior W days are all
# covered (a dietary entry, or -- with INCLUDE_ZERO_DAYS -- a 072 zero-eating day).
# This is the sample-size cost of widening the window; plotted as the blue axis to
# show the balance between cohort size and diet-microbiome fit.
window_counts <- map_dfr(Ws, function(W) {
  meta |> select(pid, sdrt, sampleid) |> tidyr::crossing(k = 1:W) |>
    mutate(fdrt = sdrt - k) |>
    left_join(diet_days, by = c("pid", "fdrt")) |>
    group_by(sampleid) |> summarise(nok = sum(!is.na(has)), .groups = "drop") |>
    summarise(W = W, n_avail = sum(nok == W))
})
message("\nPer-window available samples (blue axis):")
print(window_counts)

# ---- 6. F2a: one figure per scoring method, separate file names ------------
# Left axis: fit metric (food = solid, macro = dashed), oriented so the tightest-
# fitting window is the HIGHEST point. Right (blue) axis: number of stool samples
# available at that window (blue triangles + line).
draw_f2a <- function(df, subtitle, ylab, counts, ylim = NULL) {
  # The blue sample-count line is mapped onto an INSET of the y range so its end
  # points (esp. the day-1 maximum) are not clipped at the panel edge; the left
  # axis range is chosen tight enough that the window-to-window correlation
  # differences are visible.
  if (is.null(ylim)) {
    inset <- range(df$score)
  } else {
    sp <- diff(ylim); inset <- c(ylim[1] + 0.05 * sp, ylim[2] - 0.05 * sp)
  }
  cr <- range(counts$n_avail)
  b <- diff(inset) / diff(cr); a <- inset[1] - b * cr[1]   # map count -> y coord
  counts <- counts |> mutate(y = a + b * n_avail)
  y_scale <- if (is.null(ylim)) {
    scale_y_continuous(sec.axis = sec_axis(~ (. - a) / b, name = "n stool samples available"))
  } else {
    scale_y_continuous(limits = ylim, expand = expansion(mult = 0.02),
                       sec.axis = sec_axis(~ (. - a) / b, name = "n stool samples available"))
  }
  ggplot(df, aes(x = factor(W), y = score)) +
    geom_line(aes(group = pair, linetype = pair), linewidth = 0.3) +
    geom_point(aes(group = pair), size = 1.3) +
    geom_line(data = counts, aes(x = factor(W), y = y, group = 1),
              colour = "blue", linewidth = 0.3) +
    geom_point(data = counts, aes(x = factor(W), y = y),
               colour = "blue", shape = 17, size = 1.9) +
    y_scale +
    scale_linetype_manual(values = c("Food group based" = "solid",
                                     "Macronutrient based" = "dashed")) +
    labs(x = "diet exposure days", y = ylab, subtitle = subtitle) +
    theme_classic() +
    theme(aspect.ratio = 1.2, legend.position = "right",
          legend.title = element_blank(),
          plot.subtitle = element_text(size = 8),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          axis.title.y.right = element_text(colour = "blue"),
          axis.text.y.right  = element_text(colour = "blue"),
          axis.line.y.right  = element_line(colour = "blue"),
          axis.ticks.y.right = element_line(colour = "blue"))
}

cohort_note <- if (INCLUDE_ZERO_DAYS) "  [zero-eating days incl.]" else ""
# Symmetric: Procrustes correlation on a focused axis so the window-to-window
# differences are visible (the food curve's gain is front-loaded at 1->2 days,
# then flattens, while samples decline -- making the 2-day window the balance).
save_panel(draw_f2a(filter(procrustes_all, method == "symmetric"),
                    paste0("symmetric Procrustes; blue = n samples", cohort_note),
                    "Procrustes correlation", window_counts, ylim = c(0.15, 0.45)),
           paste0("F2a_procrustes_symmetric", suffix, ".pdf"), width = 120, height = 90)

save_panel(draw_f2a(filter(procrustes_all, method == "asymmetric"),
                    paste0("asymmetric Procrustes; blue = n samples", cohort_note),
                    "fit:  max(SS) - SS", window_counts),
           paste0("F2a_procrustes_asymmetric", suffix, ".pdf"), width = 120, height = 90)
