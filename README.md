# Reproducing the public figures of *Sugar-rich foods exacerbate antibiotic-induced microbiome injury*

This repository reproduces the **publicly-shareable** (non-PHI) figures of the dietary-sugar / microbiome manuscript from a small set of de-identified released tables. Every panel is rebuilt by a clean, numbered script under [`reproduce/`](reproduce/); the scripts read only from [`released_data/`](released_data/) and write one PDF per panel into `results/`.

The companion data deposit is on Zenodo: [**10.5281/zenodo.20278682**](https://doi.org/10.5281/zenodo.20278682) — *Supplementary Data for "Sugar-rich foods exacerbate antibiotic-induced microbiome injury"* (Dai, Anqi; CC-BY-4.0). Three of the released tables here are the same files deposited there (see [Released data](#released-data-tables)).

------------------------------------------------------------------------

## Contents

1.  [What is and isn't reproducible](#what-is-and-isnt-reproducible)
2.  [Repository layout](#repository-layout)
3.  [Environment setup](#environment-setup)
4.  [How to reproduce the figures](#how-to-reproduce-the-figures)
5.  [Scripts → figures](#scripts--figures)
6.  [Released data tables](#released-data-tables)
7.  [Running QIIME 2 through Docker](#running-qiime-2-through-docker)
8.  [Regenerating the TaxUMAP embedding](#regenerating-the-taxumap-embedding)

------------------------------------------------------------------------

## What is and isn't reproducible

**Reproducible here** (only de-identified data needed): the diet/microbiome figures — Figure 1, Figure 2, Figure 4, and Extended Data Figures E1 (c–h), E2 (b,c), E3, E4, E5, E7, E8, E9, plus the mouse experiments.

**Reproducible with restricted (internal) data**: a few panels need de-identified tables that carry **no PHI** but are **not cleared for public release**. Those tables are not shipped; they live in a gitignored `restricted_data/` folder that internal users supply. The scripts that need them are under [`reproduce/restricted/`](reproduce/restricted/) and **skip cleanly** for anyone who only has `released_data/`. Currently:

| Figure | Script | Restricted input |
|--------|--------|------------------|
| **Extended Fig. E2a** | `reproduce/restricted/60_e2a_abx_heatmap.R` | `R21_meds_updated_all_medication_classified.csv` — each patient's *full* daily antibiotic time course (the released `Data_S4` only carries the 2-day window prior to each stool sample) |

**NOT reproducible from this repository** — these depend on **patient clinical characteristics and outcomes**, which are limited-access protected health information and are *not* released:

| Figure / table | Why it cannot be rebuilt |
|--------------------------|--------------------------------------------|
| **Figure 3** | Clinical-outcome / survival analyses requiring patient outcome data |
| **Extended Fig. E1b** | "Microbiome variance explained" by clinical covariates (disease category, patient factors, transplant events) |
| **Extended Fig. E6c–j** | Clinical-outcome / survival panels |
| **Supplementary Tables 1–6** | Patient demographics / clinical characteristics tables |

Patient-level clinical variables and mortality outcomes are available via data sharing agreement per institutional policies.

------------------------------------------------------------------------

## Repository layout

```         
released_data/      de-identified input tables (the only data that ships)
restricted_data/    PHI-free but non-public inputs, internal only (gitignored, not shipped)
reproduce/          numbered scripts, one analysis "family" per number block
reproduce/restricted/  scripts that need restricted_data/ (skip cleanly when it is absent)
intermediate_data/  cached model fits, QIIME ordinations, etc. (created on first run)
results/            output PDFs, one per panel (created on first run)
```

Scripts are organised **by analysis family, not by printed figure**: an expensive model is fit and cached once, then consumed by the panels that need it.

- `10–19` human Figure 2 diversity family (+ procrustes, sensitivity, diet Faith PD)
- `20–27` Figure 1 and the Extended-Data E1/E2/E4 panels
- `30` Extended Fig. E3 (compositional PCoA)
- `40–44` taxon-abundance family (F4, E7)
- `50–54` mouse experiments (CFU, 16S, RNA-seq), driven by `00_run_mouse.R`
- `60+` (under `reproduce/restricted/`) panels requiring restricted, non-public inputs

------------------------------------------------------------------------

## Environment setup

### 1. R and packages (`renv`)

The project pins its packages with [`renv`](https://rstudio.github.io/renv/). Tested on **R 4.5.3**.

``` r
# from the project root, in R:
renv::restore()      # restores renv itself + any pinned packages

# analysis packages (install on first use, then snapshot):
renv::install(c("tidyverse","readxl","here","ggpubr","scales","ggrepel",
                "vegan","ashr","msigdbr","seriation"))
renv::install(c("bioc::Maaslin2","bioc::DESeq2","bioc::fgsea","bioc::edgeR",
                "bioc::org.Mm.eg.db","bioc::AnnotationDbi"))     # mouse seq panels
# human Bayesian + plotting family:
renv::install(c("brms","tidybayes","broom.mixed","posterior","ggtext","wesanderson",
                "cowplot","ggsci","furrr","dendextend","ggdendro","ggh4x","ggforce",
                "ggbeeswarm","lmerTest","ggpmisc"))
renv::snapshot()
```

### 2. Stan backend (for the brms models)

The Bayesian diversity / CLR models (`10`, `16`, `40`, `27`, `44`, …) use **brms** with the **cmdstanr** backend (rstan 2.32 does not compile against the current Apple clang). Install CmdStan once:

``` r
renv::install("stan-dev/cmdstanr")
cmdstanr::install_cmdstan()     # builds CmdStan (the helper expects ~/.cmdstan)
```

The fitting scripts cache every fit to `intermediate_data/*.rds` (`file =`), so they compile/sample once and are instant on rerun.

### 3. Docker (for the QIIME 2 steps)

A few scripts call QIIME 2 inside a Docker container (UniFrac / Bray-Curtis / Faith PD on the food tree). Install Docker Desktop and pull the image used by the pipeline — see [Running QIIME 2 through Docker](#running-qiime-2-through-docker).

### 4. Python / TaxUMAP (optional, Figure 1 e–h only)

The TaxUMAP embedding for F1 e–h ships precomputed (`released_data/taxumap_embedding.csv`), so you do **not** need Python to draw the figure. To regenerate it, see [Regenerating the TaxUMAP embedding](#regenerating-the-taxumap-embedding).

All public scripts resolve their input from `released_data/` (override with the `NUTRITION_DATA` env var). Scripts under `reproduce/restricted/` additionally read `restricted_data/` (override with `RESTRICTED_DATA`) and skip cleanly when it is absent.

------------------------------------------------------------------------

## How to reproduce the figures

All commands run from the project root.

### Mouse figures (F4c, E8, E9)

``` sh
Rscript reproduce/00_run_mouse.R                 # all of 50–54
RUN_SEQ=false Rscript reproduce/00_run_mouse.R   # CFU panels only, skip the slow seq fits
```

### Human figures

Run the **fitting** scripts first (they cache the models), then the **plotting** scripts that consume them. Independent families can be run in any order.

``` sh
# Figure 2 diversity family
Rscript reproduce/10_fit_diversity_models.R      # caches the F2 brms fits
Rscript reproduce/11_fig2_diversity.R            # F2 c–h, E5i
Rscript reproduce/12_diagnostics_diversity.R     # E1g, E1h
Rscript reproduce/13_fig2_procrustes.R           # F2a   (QIIME — see Docker section)

# Figure 1
Rscript reproduce/19_diet_faith_pd.R             # rebuilds diet-alpha-diversity.tsv (QIIME)
Rscript reproduce/20_fig1_taxumap.R              # F1 e–h
Rscript reproduce/21_fig1_diet_timecourse.R      # F1 b,c,i–m
Rscript -e 'rmarkdown::render("reproduce/22_fig1_food_tree.Rmd")'   # F1d (GraPhlAn/Docker)
Rscript reproduce/23_fig1_beta_diversity.R       # F1 n,o (QIIME)

# Extended Data
Rscript reproduce/24_e1cde_random_intercepts.R   # E1 c,d,e   (needs the F2d fit from 10)
Rscript reproduce/25_e1f_alpha_breakdown.R       # E1f
Rscript reproduce/26_e2bc_abx_exposure.R         # E2 b,c
Rscript reproduce/restricted/60_e2a_abx_heatmap.R  # E2a (needs restricted_data/; skips cleanly if absent)
Rscript reproduce/16_fit_e4_models.R             # caches E4 fits
Rscript reproduce/17_fig_e4.R                    # E4 a–e,i,j
Rscript reproduce/27_e4h_fndds_zscored.R         # E4h
Rscript reproduce/30_extdata_pcoa.R              # E3 a,b,c

# Taxon family (F4, E7)
Rscript reproduce/40_fit_taxon_models.R          # caches the CLR fits
Rscript reproduce/41_fig4_human.R                # F4a, F4b
Rscript reproduce/42_extdata_taxon.R             # E7 a,b,e
Rscript reproduce/43_extdata_enterococcus_asv.R  # E7 c,d
Rscript reproduce/44_e7f_abxclass.R              # E7f
```

Common environment toggles:

- `RUN_QIIME=false` — reuse the cached QIIME ordinations instead of re-running Docker.
- `QIIME2_IMAGE=...` — override the QIIME 2 image tag.
- `RUN_SEQ=false` — (mouse) skip the slow MaAsLin2 / DESeq2 / fgsea fits.

------------------------------------------------------------------------

## Scripts → figures

| Script | Panel(s) |
|----------------------------------------------------|-----------------------|
| `21_fig1_diet_timecourse.R` | **F1** b, c, i, j, k, l, m |
| `22_fig1_food_tree.Rmd` | **F1** d |
| `20_fig1_taxumap.R` | **F1** e, f, g, h |
| `23_fig1_beta_diversity.R` | **F1** n, o |
| `13_fig2_procrustes.R` | **F2** a |
| `11_fig2_diversity.R` | **F2** c, d, e, f, g, h (and E5i) |
| `41_fig4_human.R` | **F4** a, b |
| `50_figure_4c.R` | **F4** c (mouse) |
| `24_e1cde_random_intercepts.R` | **E1** c, d, e |
| `25_e1f_alpha_breakdown.R` | **E1** f |
| `12_diagnostics_diversity.R` | **E1** g, h |
| `26_e2bc_abx_exposure.R` | **E2** b, c |
| `30_extdata_pcoa.R` | **E3** a, b, c |
| `17_fig_e4.R` | **E4** a, b, c, d, e, i, j |
| `17b_e4fg_wweia.R` | **E4** f, g (WWEIA nomenclature) |
| `27_e4h_fndds_zscored.R` | **E4** h |
| `28_e5abcd_added_sugars.R` | **E5** a, b, c, d (added vs other sugars) |
| `29_e6ab_sweet_grains.R` | **E6** a, b (Sweet vs Other Grains; E6c–j are restricted) |
| `14_sensitivity_diversity.R`, `15_robustness_subsampling.R`, `18_fig_e5jk_ons.R` | **E5** e, f, g, h, j, k, l, … |
| `42_extdata_taxon.R` | **E7** a, b, e |
| `43_extdata_enterococcus_asv.R` | **E7** c, d |
| `44_e7f_abxclass.R` | **E7** f |
| `53_extdata_8_16s.R` | **E8** a–f (mouse 16S) |
| `51_extdata_8.R` | **E8** g, h, j, k, l, m (mouse) |
| `52_extdata_9.R` | **E9** a, b, c, d, h (mouse) |
| `54_extdata_9_rnaseq.R` | **E9** e, f, g (mouse RNA-seq) |

Fitting/data scripts (no panel of their own): `10`, `16`, `40` (cache the brms fits); `19` (rebuilds `diet-alpha-diversity.tsv`); `13d/13e` (procrustes sensitivity); `00_run_mouse.R` (mouse driver).

------------------------------------------------------------------------

## Released data tables

Everything in `released_data/` is de-identified and shareable. **Zenodo** column links each table to the [companion deposit](https://doi.org/10.5281/zenodo.20278682) where it is also archived; tables marked *derived* are de-identified products built in this project's upstream pipeline and are not separately on Zenodo.

### Core inputs

| File | Used for | Zenodo |
|------|----------------------------------------------------|--------|
| `152_combined_DTB.csv` | the diet tracker — every food item each patient ate; the source of all diet exposures, the food tree, and diet diversity | **same file on Zenodo** |
| `153_combined_META.csv` | the per-stool-sample analysis table (1009 samples / 158 patients): diversity outcome, antibiotic/TPN/EN exposure, prior-2-day food-group & macronutrient intake | **same file on Zenodo** |
| `Data_S4_Medication_Exposures_in_the_Two_Days_Prior_to_Stool_Sample_Collection.csv` | medications in the 2 days before each stool sample, by drug class — antibiotic-exposure panels (E2b/c) and the antibiotic-class CLR model (E7f) | **Zenodo "Data S4"** |
| `R59_meta_expanded.csv` | `153` plus the *E. faecium* CLR outcome and extra covariates — used by the taxon CLR models (F4b, E7e, E7f) | derived |
| `FPED_1516.xls`, `FPED_1720.xls` | USDA Food Patterns Equivalents Database (2015-16 and 2017-20) — the **added-sugars** content (teaspoon equivalents per 100 g) per food code, used to split total sugars into added vs other for E5a–d (`28`); 1516 is preferred, the two salad-dressing codes that exist only in 2017-20 are filled from 1720; related to Zenodo "Data S5" (FNDDS nutrient values) | reference (USDA) |
| `2015-2016 FNDDS At A Glance…xlsx`, `2019-2020 FNDDS At A Glance…xlsx` | USDA FNDDS — each food code's **WWEIA food category** (and per-100 g nutrients), used to re-derive the diversity model under WWEIA nomenclature for E4f–g (`17b`); 2015-16 preferred, 2019-20 fills the rest | Zenodo "Data S5" (FNDDS nutrient values) |

**`152_combined_DTB.csv`** — one row per food item consumed:

| column | meaning |
|------------------------------------------|------------------------------|
| `pid` | de-identified patient ID |
| `fdrt` | food day relative to transplant (day 0 = transplant) |
| `Meal`, `Food_NSC`, `Unit`, `Por_eaten` | meal, food name, unit, portions eaten |
| `Food_code`, `description` | USDA FNDDS food code and description |
| `Calories_kcal`, `Protein_g`, `Fat_g`, `Carbohydrates_g`, `Fibers_g`, `Sugars_g` | nutrient content of the portion |
| `total_weight`, `dehydrated_weight` | wet and dehydrated weight (g); `dehydrated_weight` is the intake measure used throughout |

**`153_combined_META.csv`** / **`R59_meta_expanded.csv`** — one row per stool sample:

| column | meaning |
|------------------------------------------|------------------------------|
| `pid`, `sampleid` | patient and stool-sample ID |
| `sdrt` | stool day relative to transplant |
| `simpson_reciprocal` | inverse Simpson microbiome alpha-diversity (outcome) |
| `empirical` | antibiotic-exposed in the prior 2 days (TRUE/FALSE) |
| `intensity` | conditioning intensity (nonablative / reduced / ablative) |
| `EN`, `TPN` | enteral / total parenteral nutrition exposure |
| `fg_fruit … fg_legume` | prior-2-day average intake of the 9 FNDDS food groups (g/day) |
| `ave_Calories_kcal … ave_Sugars_g` | prior-2-day average macronutrient intake |
| `timebin` | 7-day transplant-time bin (random-effect grouping) |
| *(R59 only)* `asv_1_clr`, `ci_cleaned_numeric`, `disease_lineage`, `PCA`, `exposure_type` | *E. faecium* CLR outcome + extra covariates |

**`Data_S4_Medication_Exposures…csv`** — one row per medication exposure: `sampleid`, `pid`, `sdrt`, `class`, `drug_name_clean`, `route_clean`, `drug_category_for_this_study` (broad_spectrum / fluoroquinolones / other_antibacterials / not_antibacterial).

### Microbiome tables (derived, de-identified)

| File | Columns | Used for |
|------|--------------------------|--------------------------------|
| `171_quality_asv_relab_pident97_genus.csv` | `asv_key, sampleid, count_relative, genus` | per-ASV 16S relative abundance (genus-annotated); rebuilds genus relab for F1n/o, F4a, E7b |
| `171_genus_CLR_res.csv` | `genus, sampleid, clr` | per-genus centred-log-ratio abundances for the E7a genus models |
| `171_16S_enterococcus_asv_relab.csv` | `sampleid, asv, species, relab` | 16S *Enterococcus* ASV composition for E7c (and its species-level ordering) |
| `mgx_enterococcus_species_relab.csv` | `sample, species, relab` | metagenomic *Enterococcus* species abundance for E7d |

### Diet / food-tree tables

| File | Columns / contents | Used for |
|------|--------------------------|--------------------------------|
| `072_total_patients_zero_eating_days_pid.csv` | `fdrt, diet_data_status, pid` | documented zero-intake days (distinguishes "ate nothing" from "no record") |
| `diet-alpha-diversity.tsv` | `(patient-day id), faith_pd` | diet Faith phylogenetic diversity per patient-day (regenerated by `19`) |
| `taxumap_embedding.csv` | `index_column, taxumap1, taxumap2` | precomputed TaxUMAP 2-D diet embedding for F1 e–h |
| `output_food_tree_datatree.newick` | rooted food taxonomy tree, **Food_code** leaves | UniFrac on diet (F2a, F1n, diet Faith PD) |
| `output_food_tree_datatree_name.newick` | same tree, description leaves | F1d food-tree rendering |
| `final_table_for_writing_out_to_newick.csv` | `FoodID, description, L1…L7, newickstring…, taxonomy` | the 7-level food taxonomy behind the tree (and the TaxUMAP taxonomy) |
| `NodeLabelsMCT.txt` | MCT level-code → label | food-tree node labelling |
| `food_group_color_key_final.csv` | `fgrp1, fdesc, fg1_name, color, shortname` | food-group → colour / short-name key for the forests and F1d |
| `annotation.base.txt` | GraPhlAn sector header | F1d circular-tree annotation |

### Mouse data

| File | Contents | Used for |
|------|----------------------------------------------|----------------------|
| `Dai_mouse_figure_raw_data.xlsx` | Nature source-data workbook, one sheet per mouse panel (CFU, 16S relab + metadata, RNA-seq counts) | all of F4c / E8 / E9 (scripts `50–54`) |

------------------------------------------------------------------------

## Running QIIME 2 through Docker

Scripts `13_fig2_procrustes.R`, `19_diet_faith_pd.R`, and `23_fig1_beta_diversity.R` compute phylogenetic diversity (UniFrac / Faith PD) on the **food tree** using QIIME 2. They call QIIME inside Docker so you don't need a local QIIME install.

1.  Install **Docker Desktop** and start it.

2.  Pull the image (one-time, \~10 GB):

    ``` sh
    docker pull quay.io/qiime2/qiime2:2026.4
    ```

3.  Run the script normally — it builds the feature tables in R, then runs **one** container that does `biom convert → import → diversity → export` for the whole job (one cold start, not dozens):

    ``` sh
    Rscript reproduce/13_fig2_procrustes.R
    RUN_QIIME=false Rscript reproduce/13_fig2_procrustes.R   # reuse cached ordinations
    ```

    The container is invoked as `docker run --rm --platform linux/amd64 -v <repo>:/data -w /data <image> …`, mounting the repo at `/data`.

Notes / gotchas:

- **Apple Silicon (M-series):** the QIIME image is x86-only, so it runs under `--platform linux/amd64` (Rosetta) emulation — correct but slow. That's inherent to running x86 QIIME on ARM, not the pipeline.
- **`q2-composition` / phyloseq:** the `2026.4` image's plugin manager aborts trying to import the R package `phyloseq` (not installed) via `q2-composition`. The scripts run `pip uninstall -y q2-composition` at container start (the container is `--rm`, so it's ephemeral and `q2-composition` is unused here).
- **Figure 1d (`22_fig1_food_tree.Rmd`)** uses a *different* container — **GraPhlAn** via the `shengwei/graphlan` image (GraPhlAn runs only under Python 2.7) — for the circular food-taxonomy render. The R steps build the inputs; the notebook's Docker commands draw the panel.

------------------------------------------------------------------------

## Regenerating the TaxUMAP embedding

Figure 1 e–h is drawn by `20_fig1_taxumap.R` from the precomputed `released_data/taxumap_embedding.csv`, so **the figure reproduces without Python**. To regenerate the embedding from scratch, the full step-by-step (build the two input tables in R → install TaxUMAP → run it) is in [**`reproduce/taxumap_pipeline_HOWTO.md`**](reproduce/taxumap_pipeline_HOWTO.md).

In brief:

``` sh
# 1. install TaxUMAP (https://github.com/jsevo/taxumap) in its own env
git clone https://github.com/jsevo/taxumap.git && cd taxumap
conda create -n taxumap python=3.9 -y && conda activate taxumap
pip install -e .

# 2. (in R) build food_code_relative.csv + food_taxa.csv from 152_combined_DTB.csv
#    and final_table_for_writing_out_to_newick.csv  — see the HOWTO

# 3. run TaxUMAP with the manuscript's weights
run_taxumap.py \
  -t intermediate_data/food_taxa.csv \
  -m intermediate_data/food_code_relative.csv \
  --agg_levels Kingdom/Phylum/Class/Order \
  -n 173 --weights 0.01/2/10/10 \
  -o intermediate_data
```

UMAP is stochastic (the embedding can rotate/flip between runs), so for an **exact** figure match use the shipped `taxumap_embedding.csv` rather than a fresh embedding.

------------------------------------------------------------------------

*Patient-level clinical variables and mortality outcomes (Figure 3, E1b, E6c–j, Supplementary Tables 1–6) are available via data sharing agreement per institutional policies. (E2a is reproducible internally from the restricted antibiotic time-course table, see [What is and isn't reproducible](#what-is-and-isnt-reproducible).)*
