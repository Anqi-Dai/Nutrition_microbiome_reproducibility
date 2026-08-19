# Reproducing the public figures of *Sugar-rich foods exacerbate antibiotic-induced microbiome injury*

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21290618-blue.svg)](https://doi.org/10.5281/zenodo.21290618) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository reproduces the **publicly-shareable** figures — those free of protected health information (PHI) — of the dietary-sugar / microbiome manuscript from a small set of de-identified released tables. Every panel is rebuilt by a clean, numbered script under [`reproduce/`](reproduce/); the scripts read only from [`released_data/`](released_data/) and write one PDF per panel into `results/`.

The companion data deposit is on Zenodo: [**10.5281/zenodo.14538105**](https://doi.org/10.5281/zenodo.14538105) — *Supplementary Data for "Sugar-rich foods exacerbate antibiotic-induced microbiome injury"* (Dai, Anqi; CC-BY-4.0). Three of the released tables here are the same files deposited there (see [Released data](#released-data-tables)).

------------------------------------------------------------------------

## Contents

1.  [All data needed to reproduce the figures is released publicly](#all-data-needed-to-reproduce-the-figures-is-released-publicly)
2.  [Repository layout](#repository-layout)
3.  [Environment setup](#environment-setup)
4.  [A note on the former restricted tier](#a-note-on-the-former-restricted-tier)
5.  [How to reproduce the figures](#how-to-reproduce-the-figures)
6.  [Scripts → figures](#scripts--figures)
7.  [Released data tables](#released-data-tables)
8.  [Running QIIME 2 through Docker](#running-qiime-2-through-docker)
9.  [Regenerating the TaxUMAP embedding](#regenerating-the-taxumap-embedding)
10. [Citation](#citation)
11. [License](#license)
12. [Acknowledgements](#acknowledgements)

------------------------------------------------------------------------

## All data needed to reproduce the figures is released publicly

Every panel in the paper is reproducible from the de-identified tables shipped in
[`released_data/`](released_data/). Clone the repo, restore the environment, and run —
there is no second tier to request and no access step.

This changed on 2026-08-19. The clinical-outcome and antibiotic-exposure tables behind
Fig. 3a,b and Extended Data Figs 1b, 2a and 6c–j were previously held back as
"PHI-free but not cleared for release" and fetched separately by internal users. They
have since been cleared, reduced to the columns the analyses actually use, and moved
into `released_data/`. Nothing needs to be requested from the authors any more; if you
have read an older version of this README or the data availability statement, that
instruction is out of date.

The two newly released tables:

| Table | What it is |
|-------|------------|
| `df_main_clinical_outcome.csv` | 173 patients, one row each, 20 columns: transplant characteristics, day-12-landmarked overall survival, engraftment and discharge landmarks, broad-spectrum antibiotic exposure, diet-pattern cluster, and the peri-transplant diet summary |
| `R21_meds_updated_all_medication_classified.csv` | 95,248 medication courses, one row each: drug, route, start and stop day relative to transplant, and the per-study antibiotic class. Each patient's *complete* daily record — a superset of the released `Data_S4`, which carries only the two-day window before each stool sample |

Both are de-identified: patients appear only as the arbitrary `pid` codes already used
throughout the release, there are no dates, and all times are days relative to
transplant.

The scripts that read them are numbered 60–67 in [`reproduce/`](reproduce/),
alongside every other panel script.

| Figure | Script |
|--------|--------|
| **Extended Fig. E2a** | `reproduce/60_e2a_abx_heatmap.R` — per-patient daily antibiotic-exposure heatmap |
| **Data S6** | `reproduce/61_dataS6_pt_timecourse.R` — per-patient calorie / diet-diversity / microbiome-diversity timecourse; the clinical table supplies the engraftment day |
| **Fig. 3 a,b · E6 c,d,j · Supp. Tables 1–7** | `reproduce/63_fig3_e6_clinical.R` — the survival models, the sugar-vs-calorie scatters, the HR contour, and the supplementary characteristics tables. Supplementary Table 7 (mouse diet composition) is built by `reproduce/51_extdata_8.R` and read from its cache in `intermediate_data/`, so run 51 first to get `Supplementary_Tables_1_7.pdf`; without it you get `Supplementary_Tables_1_6.pdf` |
| **Extended Fig. E6 e,f,g** | `reproduce/64_e6efg_cluster_intake.R` — daily calorie / macronutrient intake by diet-pattern cluster |
| **Extended Fig. E6h** | `reproduce/66_e6h_alpha_trajectory.R` — fecal alpha-diversity trajectory by cluster |
| **Extended Fig. E6i** | `reproduce/65_e6i_discharge.R` — cumulative incidence of hospital discharge after engraftment by cluster (adjusted HR 1.54, p = 0.023) |
| **Extended Fig. E1b** | `reproduce/67_e1b_covariates_contribution.R` — per-covariate microbiome variance explained (`vegan::envfit` r²) |

------------------------------------------------------------------------

## Repository layout

```         
released_data/      de-identified input tables (every input the figures need)
reproduce/          numbered scripts, one analysis "family" per number block
intermediate_data/  cached model fits, QIIME ordinations, etc. (created on first run)
results/            output PDFs, one per panel (created on first run)
```

Scripts are organised **by analysis family, not by printed figure**: an expensive model is fit and cached once, then consumed by the panels that need it.

- `10–19` human Figure 2 diversity family (+ procrustes, sensitivity, diet Faith PD)
- `20–27` Figure 1 and the Extended-Data E1/E2/E4 panels
- `30` Extended Fig. E3 (compositional PCoA)
- `40–44` taxon-abundance family (F4, E7)
- `50–54` mouse experiments (CFU, 16S, RNA-seq), driven by `00_run_mouse.R`
- `60–67` clinical-outcome and antibiotic-exposure panels (Fig. 3, E1b, E2a, E6 c–j, Data S6, Supp. Tables 1–7)

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

The TaxUMAP embedding for F1 e–h ships precomputed (`released_data/taxumap_embedding.csv`), so you do not need Python to draw the figure. To regenerate it, see [Regenerating the TaxUMAP embedding](#regenerating-the-taxumap-embedding).

Every script resolves its input from `released_data/` (override with the `NUTRITION_DATA` env var), including those under `reproduce/`.

------------------------------------------------------------------------

## A note on the former restricted tier

Earlier versions of this repository kept two tables outside `released_data/`, in a
separate non-public tier that lab members fetched from an internal share and everyone
else had to request under a data-sharing agreement. Some panel scripts lived in their own
`reproduce/restricted/` folder and skipped when those inputs were absent.

**None of that applies any more.** Both tables were cleared for public release on
2026-08-19 and ship in `released_data/`; the scripts moved into `reproduce/` with
everything else; and there is no separate fetch step, no access request, and no skipping.
A clean clone reproduces every panel. If you are following older instructions that
mention a restricted tier, they are out of date.

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
Rscript reproduce/26b_e2de_fluoroquinolone.R     # E2 d,e   (fits the simplified fluoroquinolone model)
Rscript reproduce/60_e2a_abx_heatmap.R  # E2a (per-patient daily antibiotic exposure)
Rscript reproduce/61_dataS6_pt_timecourse.R  # Data S6 (per-patient timecourse)
Rscript reproduce/63_fig3_e6_clinical.R      # Fig 3 a,b + E6 c,d,j + Supp. Tables 1-7 as one PDF (reads the cleaned df_main; S7 comes from 51's cache)
Rscript reproduce/64_e6efg_cluster_intake.R  # E6 e,f,g (daily intake by diet-pattern cluster)
Rscript reproduce/65_e6i_discharge.R         # E6i (discharge cumulative incidence by cluster, adjusted HR)
Rscript reproduce/66_e6h_alpha_trajectory.R  # E6h (fecal alpha-diversity trajectory by cluster)
Rscript reproduce/67_e1b_covariates_contribution.R  # E1b (per-covariate microbiome variance explained)
Rscript reproduce/16_fit_e4_models.R             # caches E4 fits
Rscript reproduce/17_fig_e4.R                    # E4 a–e,i,j
Rscript reproduce/27_e4h_fndds_zscored.R         # E4h
Rscript reproduce/30_extdata_pcoa.R              # E3 a,b,c

# Taxon family (F4, E7)
Rscript reproduce/45_quality_genus_relab.R       # rebuilds 45_quality_asv_relab_pident97_genus.csv from the 63 tables
Rscript reproduce/40_fit_taxon_models.R          # caches the CLR fits (derives genus CLR inline)
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
| `26b_e2de_fluoroquinolone.R` | **E2** d, e (simplified diversity model with an explicit prophylactic-fluoroquinolone term) |
| `30_extdata_pcoa.R` | **E3** a, b, c |
| `17_fig_e4.R` | **E4** a, b, c, d, e, i, j |
| `17b_e4fg_wweia.R` | **E4** f, g (What We Eat in America, WWEIA, nomenclature) |
| `27_e4h_fndds_zscored.R` | **E4** h |
| `28_e5abcd_added_sugars.R` | **E5** a, b, c, d (added vs other sugars) |
| `29_e6ab_sweet_grains.R` | **E6** a, b (Sweet vs Other Grains; E6c–j come from 63/64/65/66) |
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

Everything in `released_data/` is de-identified and shareable. **Zenodo** column links each table to the [companion deposit](https://doi.org/10.5281/zenodo.14538105) where it is also archived; tables marked *derived* are de-identified products built in this project's upstream pipeline and are not separately on Zenodo.

### Core inputs

| File | Used for | Zenodo |
|------|----------------------------------------------------|--------|
| `152_combined_DTB.csv` | the diet tracker — every food item each patient ate; the source of all diet exposures, the food tree, and diet diversity | **same file on Zenodo** |
| `153_combined_META.csv` | the per-stool-sample analysis table (1009 samples / 158 patients): diversity outcome, antibiotic/TPN/EN exposure, prior-2-day food-group & macronutrient intake | **same file on Zenodo** |
| `Data_S4_Medication_Exposures…csv` | medications in the 2 days before each stool sample, by drug class — antibiotic-exposure panels (E2b/c), the fluoroquinolone diversity model (E2d/e) and the antibiotic-class CLR model (E7f) | **Zenodo "Data S4"** |
| `R59_meta_expanded.csv` | `153` plus the *E. faecium* CLR outcome and extra covariates — used by the taxon CLR models (F4b, E7e, E7f) | derived |
| `FPED_1516.xls`, `FPED_1720.xls` | USDA Food Patterns Equivalents Database (2015-16 and 2017-20) — the **added-sugars** content (teaspoon equivalents per 100 g) per food code, used to split total sugars into added vs other for E5a–d (`28`); 1516 is preferred, the two salad-dressing codes that exist only in 2017-20 are filled from 1720; related to Zenodo "Data S5", the nutrient values from the USDA Food and Nutrient Database for Dietary Studies (FNDDS) | reference (USDA) |
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

**`df_main_clinical_outcome.csv`** — one row per patient (173), the clinical-outcome table behind Fig. 3a,b and Extended Data Figs 1b, 6c–j. Released 2026-08-19; reduced to the 20 columns the analyses read:

| column | meaning |
|------------------------------------------|------------------------------|
| `pid` | de-identified patient ID (same codes as the tables above) |
| `age`, `sex`, `disease.simple` | age at transplant (years), sex, simplified primary disease |
| `source` | graft: Unmodified / TCD (T-cell-depleted) / Cord |
| `intensity` | conditioning intensity: nonablative / reduced / ablative |
| `gvhd_ppx`, `source_and_gvhdppx` | GVHD prophylaxis, and its combination with graft source (the survival-model adjustment) |
| `OStime_30`, `OSevent` | overall survival landmarked at day 12, and the death indicator. `OStime_30` is a misnomer for `OStime_12`, kept to match the analysis code |
| `engraftment_day` | neutrophil engraftment day (ANC 500); NA for 5 patients |
| `tLOS`, `LOS` | days from engraftment to discharge, and whether discharge was observed (the E6i outcome) |
| `leng_of_stay` | total hospital length of stay — a different quantity from `tLOS`, and the one in Supp. Tables S1/S2 |
| `modal_diet` | diet-pattern cluster: Cluster 1 (n=114) / Cluster 2 (n=59) |
| `day_exposed` | days on any broad-spectrum antibiotic within HCT day −7…+12 (0–20) |
| `SugarCal_cat_high` | above/below cohort-median sugar density |
| `avg_total_caloric_intake`, `avg_sugar_daily_intake_gram`, `avg_sugar_density_per_1000kcal` | means per day over HCT day −7…+12: kcal, g sugar, g sugar per 1000 kcal |

Five columns are categorical and **their level order is fixed, not alphabetical**, because the first level is the reference category in every model built on this table. Read the file with `read_clinical()` in [`reproduce/human/_human_helpers.R`](reproduce/human/_human_helpers.R), which restores them and stops if a category is unexpected or missing; if you re-analyse the table in your own code, set the levels yourself.

**`R21_meds_updated_all_medication_classified.csv`** — one row per medication course (95,248), each patient's *complete* daily record and the source of `day_exposed` above. Released 2026-08-19. Columns: `pid`, `drug_name_clean` (382 distinct drugs), `drug_route` as `{drug}__{route}`, `startday` and `stopday` in days relative to transplant (inclusive), and `drug_category_for_this_study` (broad_spectrum / fluoroquinolones / other_antibacterials / not_antibacterial). This is a superset of `Data_S4` above, which carries only the two-day window before each stool sample.

### Microbiome tables (derived, de-identified)

| File | Columns | Used for |
|------|--------------------------|--------------------------------|
| `63_asv_count_relab_res.csv` | `asv_key, sampleid, count, count_relative` | per-ASV 16S counts and whole-community relative abundance (one row per observed, non-zero count); the ASV-count source for E1b, and joined to the annotation for E7c and the quality-genus / genus-CLR tables |
| `63_asv_blast_annotation.csv` | `asv_key` + taxonomy + BLAST hit (see below) | per-ASV taxonomic lineage and the BLAST hit that assigned it; one row per ASV |
| `45_quality_asv_relab_pident97_genus.csv` | `asv_key, sampleid, count_relative, genus` | per-ASV 16S relative abundance with a *quality* genus (kept only where BLAST `pident > 97`, else NA); rebuilds genus relab for F1n/o, F4a, E7b, the E3 PCoA and the E7a genus models |
| `mgx_enterococcus_species_relab.csv` | `sample, species, relab` | metagenomic *Enterococcus* species abundance for E7d |

**`63_asv_blast_annotation.csv` columns.** Each ASV sequence is the *query*, aligned against a 16S reference database with BLAST; the row records the assigning hit. Columns marked *BLAST* are standard BLAST tabular (outfmt 6) fields; the rest are derived by the annotation pipeline.

| Column | Meaning |
|--------|---------|
| `asv_key` | ASV identifier (e.g. `asv_1`) |
| `kingdom … species` | assigned taxonomic lineage (`kingdom, phylum, class, order, family, genus, species`); `order` was the raw export's `ordr` |
| `query_length` | length (bp) of the ASV query sequence *(BLAST `qlen`)* |
| `align_length` | length of the query–subject alignment, gaps included *(BLAST `length`)* |
| `pident` | percent of identical bases over the alignment *(BLAST `pident`)* |
| `nident` | number of identical bases in the alignment *(BLAST `nident`)* |
| `score` | raw alignment score *(BLAST `score`)* |
| `length_ratio` | `align_length / query_length`: fraction of the query covered by the alignment (~1 = full-length match) *(derived)* |

### Diet / food-tree tables

| File | Columns / contents | Used for |
|------|--------------------------|--------------------------------|
| `072_total_patients_zero_eating_days_pid.csv` | `fdrt, diet_data_status, pid` | documented zero-intake days (distinguishes "ate nothing" from "no record") |
| `diet-alpha-diversity.tsv` | `(patient-day id), faith_pd` | diet Faith phylogenetic diversity per patient-day (regenerated by `19`) |
| `taxumap_embedding.csv` | `index_column, taxumap1, taxumap2` | precomputed TaxUMAP 2-D diet embedding for F1 e–h |
| `output_food_tree_datatree.newick` | rooted food taxonomy tree, **Food_code** leaves | UniFrac on diet (F2a, F1n, diet Faith PD) |
| `output_food_tree_datatree_name.newick` | same tree, description leaves | F1d food-tree rendering |
| `final_table_for_writing_out_to_newick.csv` | `FoodID, description, L1…L7, newickstring…, taxonomy` | the 7-level food taxonomy behind the tree (and the TaxUMAP taxonomy) |
| `NodeLabelsMCT.txt` | `Level.code, Main.food.description` (the label), plus `n_digits, source, usda_fndds_description, used_in_F1d_tree` (provenance) | food-tree node labelling. Levels 1–3 are official USDA FNDDS subgroups; levels 4–6 are mostly the finer split from the Knights lab food-tree code, which USDA does not publish. See the provenance section of `22_fig1_food_tree.Rmd` for the citations |
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

## Citation

Machine-readable metadata lives in [`CITATION.cff`](CITATION.cff) — use GitHub's **Cite this repository** button.

The accompanying manuscript is **under review**. Once it is published, please **cite the article in preference to this software**; this repository will then carry a `preferred-citation` entry pointing at it.

Two related DOIs, which are *not* interchangeable:

| DOI | What it is |
|-----------------------------|------------------------------------------|
| [10.5281/zenodo.21290618](https://doi.org/10.5281/zenodo.21290618) | **Software (all versions)** — the *concept* DOI. It always resolves to the latest release. **Cite this one.** |
| [10.5281/zenodo.14538105](https://doi.org/10.5281/zenodo.14538105) | **Data (all versions)** — the companion Supplementary Data deposit's *concept* DOI. **Cite this one.** |

------------------------------------------------------------------------

## License

- **Code** (everything under `reproduce/`, and the repository as a whole): [MIT](LICENSE).
- **Data** (`released_data/`): [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/), matching the [companion Zenodo deposit](https://doi.org/10.5281/zenodo.14538105).

Third-party reference tables redistributed in `released_data/` (the USDA FNDDS "At A Glance" and FPED releases) are US Government works in the public domain.

------------------------------------------------------------------------

## Acknowledgements

This repository reproduces analyses developed with several colleagues; the clean scripts here are ports of their original work:

- **William Jogia** (Institute for Systems Genetics and Department of Microbiology, NYU Langone Health Grossman School of Medicine) — the WWEIA-nomenclature food-group analysis (E4 f,g), the added-sugars analysis (E5 a–d), and the sweet-grain split (E6 a,b).
- **Mirae Baichoo** (Adult Bone Marrow Transplantation Service, Department of Medicine, Memorial Sloan Kettering Cancer Center) — the per-covariate microbiome variance-explained analysis (E1b).
- **Teng Fei** (Department of Epidemiology and Biostatistics, Memorial Sloan Kettering Cancer Center) — the clinical-outcome and survival analyses (Figure 3, and E6 h,i,j).
- **Nicholas R. Waters** (Adult Bone Marrow Transplantation Service, Department of Medicine, Memorial Sloan Kettering Cancer Center) — the *Enterococcus* ASV / species analyses establishing that ASV 1 is *E. faecium* (E7 c,d).

------------------------------------------------------------------------

*Patient-level clinical variables and mortality outcomes — underlying Figure 3, E1b, E2a, E6 c–j and Supplementary Tables 1–6 — were released publicly on 2026-08-19 and ship in `released_data/`; no data sharing agreement is required.*
