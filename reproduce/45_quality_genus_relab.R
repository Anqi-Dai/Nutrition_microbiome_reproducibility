# Build released_data/45_quality_asv_relab_pident97_genus.csv from the two 63 ASV
# tables. Per-ASV whole-community relative abundance (count_relative) carried with its
# genus, but the genus is kept only where the BLAST identity supports it (pident > 97);
# below that the genus is set to NA ("not quality"). This genus-annotated ASV relab is
# what the diet/microbiome procrustes (13), beta-diversity (23), E3 PCoA (30), the E7a
# genus models (40) and the F4a/E7b Spearman helper read.
#
# Derived exactly (to the row) from:
#   63_asv_count_relab_res.csv    asv_key, sampleid, count, count_relative
#   63_asv_blast_annotation.csv   asv_key, ..., genus, ..., pident, ...
# and reproduces the previously shipped 171_quality_asv_relab_pident97_genus.csv.
#
# Run from the repo root. Writes to the NUTRITION_DATA root (default released_data/),
# after comparing against the previously shipped file and printing whether it
# reproduces it.

suppressPackageStartupMessages(library(tidyverse))

data_root <- Sys.getenv("NUTRITION_DATA", unset = here::here("released_data"))
released <- function(f) file.path(data_root, f)

counts <- read_csv(released("63_asv_count_relab_res.csv"), show_col_types = FALSE) |>
  select(asv_key, sampleid, count_relative)

# pident is from the detailed BLAST hit; ASVs with no hit (pident NA) are not quality,
# so their genus is dropped as well.
annotation <- read_csv(released("63_asv_blast_annotation.csv"), show_col_types = FALSE) |>
  select(asv_key, genus, pident)

quality_genus <- counts |>
  left_join(annotation, by = "asv_key") |>
  mutate(genus = if_else(!is.na(pident) & pident > 97, genus, NA_character_)) |>
  select(asv_key, sampleid, count_relative, genus)

message(nrow(quality_genus), " rows; ", sum(!is.na(quality_genus$genus)),
        " with a quality (pident > 97) genus, over ",
        n_distinct(quality_genus$sampleid), " samples")

# Reproduction check against the previously shipped table when present.
old_file <- released("171_quality_asv_relab_pident97_genus.csv")
if (file.exists(old_file)) {
  old <- read_csv(old_file, show_col_types = FALSE)
  cmp <- quality_genus |>
    full_join(old, by = c("asv_key", "sampleid"), suffix = c("_new", "_old"))
  genus_mismatch <- sum(!(cmp$genus_new == cmp$genus_old |
                          (is.na(cmp$genus_new) & is.na(cmp$genus_old))), na.rm = TRUE)
  message("REPRODUCTION CHECK vs 171_quality: ",
          sum(!is.na(cmp$count_relative_new) & !is.na(cmp$count_relative_old)), "/", nrow(old),
          " rows matched, max relab diff ",
          signif(max(abs(cmp$count_relative_new - cmp$count_relative_old), na.rm = TRUE), 3),
          ", genus mismatches ", genus_mismatch)
}

out_file <- released("45_quality_asv_relab_pident97_genus.csv")
write_csv(quality_genus, out_file)
message("Wrote ", out_file)
