# Regenerating the TaxUMAP embedding (Figure 1 e-h)

`reproduce/20_fig1_taxumap.R` plots panels F1 e-h from a precomputed TaxUMAP
embedding, `released_data/taxumap_embedding.csv`. That embedding is produced once
by the Python tool **TaxUMAP** (https://github.com/jsevo/taxumap); it is shipped so
the figure reproduces without anyone installing the Python stack. This document is
for regenerating it from scratch — the upstream (Stages 1-3) half of the original
`162_taxUMAP__Figure1_e_to_h.Rmd`.

The flow is: build two input tables in R, install TaxUMAP, run it on those tables.

```
152_combined_DTB.csv ─┐
                      ├─(R, Stage 1)→ food_code_relative.csv ─┐
final_table_..._newick.csv ─(R, Stage 2)→ food_taxa.csv ──────┼─(TaxUMAP)→ taxumap_embedding.csv
                                                              ┘
```

TaxUMAP treats each patient-day as a "sample", each food code as an "ASV/OTU", and
its daily relative weight as the "abundance"; the food taxonomy (Kingdom/Phylum/
Class/Order from the food tree) plays the role of microbial taxonomy.

## Stage 1 — daily food-code relative-abundance matrix

TaxUMAP's microbiota table needs one row per sample (`index_column`) and one column
per food code holding compositional (relative) abundances. The diet tracker records
absolute dehydrated weight per item, so for each `(patient, day)` sum to a daily
total per food code, then divide by the daily grand total.

```r
library(tidyverse)
dtb <- read_csv("released_data/152_combined_DTB.csv")

fc_table <- dtb |>
  select(pid, fdrt, Food_code, dehydrated_weight) |>
  group_by(pid, fdrt, Food_code) |>
  summarise(total = sum(dehydrated_weight), .groups = "drop") |>
  mutate(index_column = str_glue("P{pid}_{fdrt}")) |>
  select(index_column, Food_code, total)

fc_total <- fc_table |>
  group_by(index_column) |>
  summarise(daily_total = sum(total), .groups = "drop")

fc_df <- fc_table |>
  left_join(fc_total, by = "index_column") |>
  mutate(daily_relative = total / daily_total) |>
  select(index_column, Food_code, daily_relative) |>
  pivot_wider(names_from = Food_code, values_from = daily_relative, values_fill = 0)

stopifnot(all(abs(rowSums(fc_df[ , -1]) - 1) < 1e-8))   # each row sums to 1
write_csv(fc_df, "intermediate_data/food_code_relative.csv")
```

## Stage 2 — food taxonomy table

TaxUMAP's taxonomy table is indexed by the ASV/OTU labels (here food codes) with
columns of higher taxonomic groups ordered left-to-right by **decreasing** hierarchy
(Kingdom -> Order). Unknown levels must be the string `nan`. Restrict to the food
codes that actually appear in the abundance matrix.

```r
actual_foodids <- setdiff(colnames(fc_df), "index_column")

food_taxa <- read_csv("released_data/final_table_for_writing_out_to_newick.csv") |>
  mutate(FoodID = as.character(FoodID)) |>
  filter(FoodID %in% actual_foodids) |>
  select(ASV = FoodID, Kingdom = L1, Phylum = L2, Class = L3, Order = L4)

# any food code without taxonomy is dropped by TaxUMAP
setdiff(actual_foodids, food_taxa$ASV)
write_csv(food_taxa, "intermediate_data/food_taxa.csv")
```

## Install TaxUMAP

TaxUMAP is a Python CLI (no PyPI/Bioconda release yet), so install it from source in
its own environment. From the README: `git clone` then `pip install -e .`.

```sh
git clone https://github.com/jsevo/taxumap.git
cd taxumap

# isolated environment (conda shown; a venv works too)
conda create -n taxumap python=3.9 -y
conda activate taxumap

# editable/developer install; pulls umap-learn, pandas, numpy, etc. and puts
# run_taxumap.py on the PATH
pip install -e .
```

Verify the CLI is available:

```sh
run_taxumap.py --help
```

## Stage 3 — run TaxUMAP

With the env active, run it on the two tables from Stages 1-2. This is the exact
command (and weights) used for the manuscript figure:

```sh
run_taxumap.py \
  -t intermediate_data/food_taxa.csv \
  -m intermediate_data/food_code_relative.csv \
  --agg_levels Kingdom/Phylum/Class/Order \
  -n 173 \
  --weights 0.01/2/10/10 \
  -o intermediate_data
```

Arguments:

| flag | value | meaning |
|------|-------|---------|
| `-t` / `--taxonomy`   | `food_taxa.csv`          | taxonomy table (Stage 2) |
| `-m` / `--microbiota` | `food_code_relative.csv` | abundance matrix (Stage 1) |
| `--agg_levels`        | `Kingdom/Phylum/Class/Order` | taxonomic levels aggregated and embedded |
| `-n` / `--neigh`      | `173`                    | UMAP `n_neighbors`; the README recommends matching the unique patient count (173 here) |
| `--weights`           | `0.01/2/10/10`           | per-level weights (one per `--agg_levels` entry), up-weighting the coarser food groups |
| `-o` / `--outdir`     | `intermediate_data`      | output directory |

TaxUMAP writes the 2-D embedding (one row per patient-day with `index_column`,
`taxumap1`, `taxumap2`). Rename/copy it to `released_data/taxumap_embedding.csv`,
and `reproduce/20_fig1_taxumap.R` will draw F1 e-h from it.

Notes:
- UMAP is stochastic; set the RNG seed (or accept that the embedding rotates/flips
  between runs). The shipped `taxumap_embedding.csv` is the manuscript's run, so for
  an exact figure match use it rather than a fresh embedding.
- `-n 173` matches the cohort's patient count; change it if you rerun on a different
  cohort.
