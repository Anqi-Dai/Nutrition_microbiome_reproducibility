# Reproduce released_data/diet-alpha-diversity.tsv: per patient-day diet alpha
# diversity (Faith's phylogenetic diversity) on the food tree.
#
# Pipeline (refactor of 005_diet_alpha_div.Rmd + the food_tree Snakefile / the
# qiime steps in 156_clinical_data.Rmd / 05_food_diversity.Rmd):
#   1. food-code count matrix from the diet tracker: rows = Food_code, columns =
#      one per patient-day (fid = "{pid}d{fdrt}", e.g. P100d-1), values = daily
#      summed dehydrated weight.
#   2. QIIME2 (one Docker container): tsv -> biom -> qza; import the rooted food
#      tree; `qiime diversity alpha-phylogenetic --p-metric faith_pd`; export.
#   3. the exported alpha-diversity.tsv IS diet-alpha-diversity.tsv (id, faith_pd).
#
# The script writes results to released_data/diet-alpha-diversity.tsv, but first
# compares the freshly computed values against the shipped file and prints whether
# they reproduce it exactly (the existing file is read into memory before any
# overwrite, so the comparison is against the original).
#
# RUN_QIIME=false reuses the cached export; QIIME2_IMAGE picks the tag.

suppressPackageStartupMessages(library(tidyverse))

ROOT      <- here::here()
released  <- function(f) file.path(ROOT, "released_data", f)
WORK      <- file.path(ROOT, "intermediate_data", "diet_faith_pd")
TREE_NWK  <- released("output_food_tree_datatree.newick")     # Food_code leaves, rooted
TARGET    <- released("diet-alpha-diversity.tsv")
QIIME2_IMAGE <- Sys.getenv("QIIME2_IMAGE", unset = "quay.io/qiime2/qiime2:2026.4")
RUN_QIIME    <- tolower(Sys.getenv("RUN_QIIME", unset = "true")) %in% c("true", "1", "yes")

dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(TREE_NWK))

# ---- 1. food-code count matrix (rows = Food_code, cols = patient-day) -------
dtb <- read_csv(released("152_combined_DTB.csv"), show_col_types = FALSE)

fcts <- dtb |>
  mutate(fid = str_glue("{pid}d{fdrt}")) |>
  group_by(fid, Food_code) |>
  summarise(daily_de_wt = sum(dehydrated_weight), .groups = "drop") |>
  pivot_wider(names_from = fid, values_from = daily_de_wt, values_fill = 0) |>
  rename(`#OTU ID` = Food_code)

matrix_tsv <- file.path(WORK, "005_food_code_counts_matrix.tsv")
write_tsv(fcts, matrix_tsv)
message(sprintf("food-code matrix: %d food codes x %d patient-days",
                nrow(fcts), ncol(fcts) - 1))

# ---- 2. QIIME2 Faith PD in one container -----------------------------------
run_qiime <- function() {
  rel <- function(p) sub(paste0("^", ROOT, "/?"), "", p)        # host path -> /data path
  base_c <- file.path("/data", rel(WORK))
  tree_c <- file.path("/data", rel(TREE_NWK))

  script <- sprintf('
set -euo pipefail
# 2026.4 image: q2-composition aborts the plugin manager importing R phyloseq; drop
# it (unused; container is --rm so this is ephemeral).
pip uninstall -y q2-composition >/dev/null 2>&1 || true
BASE=%s
biom convert -i "$BASE/005_food_code_counts_matrix.tsv" -o "$BASE/food_counts.biom" --to-hdf5 --table-type="Table"
qiime tools import --input-path "$BASE/food_counts.biom" --output-path "$BASE/food_counts.qza" --type "FeatureTable[Frequency]"
[ -f "$BASE/food_tree.qza" ] || qiime tools import --input-path %s --output-path "$BASE/food_tree.qza" --type "Phylogeny[Rooted]"
qiime diversity alpha-phylogenetic --i-table "$BASE/food_counts.qza" --i-phylogeny "$BASE/food_tree.qza" --p-metric faith_pd --o-alpha-diversity "$BASE/faith_pd.qza"
rm -rf "$BASE/faith_export"
qiime tools export --input-path "$BASE/faith_pd.qza" --output-path "$BASE/faith_export"
', base_c, tree_c)

  cmd <- sprintf('docker run --rm --platform linux/amd64 -v %s:/data -w /data %s bash -c %s',
                 shQuote(ROOT), shQuote(QIIME2_IMAGE), shQuote(script))
  message("running QIIME2 alpha-phylogenetic (faith_pd) in one container ...")
  if (system(cmd) != 0) stop("QIIME2 container exited non-zero")
}

export_tsv <- file.path(WORK, "faith_export", "alpha-diversity.tsv")
if (RUN_QIIME || !file.exists(export_tsv)) run_qiime() else message("reusing cached export")
stopifnot(file.exists(export_tsv))

# ---- 3. compare to the shipped file, then write ----------------------------
new <- read_tsv(export_tsv, show_col_types = FALSE)
names(new)[1] <- "id"

if (file.exists(TARGET)) {
  old <- read_tsv(TARGET, show_col_types = FALSE); names(old)[1] <- "id"
  cmp <- inner_join(old, new, by = "id", suffix = c("_old", "_new"))
  max_abs  <- max(abs(cmp$faith_pd_old - cmp$faith_pd_new))
  max_rel  <- max(abs(cmp$faith_pd_old - cmp$faith_pd_new) /
                    pmax(abs(cmp$faith_pd_old), 1e-9))
  byte_id  <- identical(readLines(export_tsv), readLines(TARGET))
  message("\n================  REPRODUCTION CHECK (diet-alpha-diversity.tsv)  ================")
  message(sprintf("shipped rows: %d   reproduced rows: %d   ids in common: %d",
                  nrow(old), nrow(new), nrow(cmp)))
  message(sprintf("ids only in shipped: %d   ids only in reproduced: %d",
                  length(setdiff(old$id, new$id)), length(setdiff(new$id, old$id))))
  message(sprintf("max abs diff in faith_pd: %.3e   max rel diff: %.3e", max_abs, max_rel))
  message(sprintf("byte-identical to shipped file: %s", byte_id))
  if (byte_id) {
    message(">>> EXACT MATCH: reproduction is byte-identical to the shipped file.")
  } else if (max_abs < 1e-6) {
    message(">>> NUMERIC MATCH: values agree to < 1e-6 (formatting/row-order only).")
  } else {
    message(">>> DIFFERENCE: values differ -- FLAG FOR ANGEL TO CHECK (see stats above).")
  }
  message("===============================================================================\n")
}

file.copy(export_tsv, TARGET, overwrite = TRUE)
message("wrote ", TARGET)
