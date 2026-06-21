# Mouse figure reproduction roadmap

Covers every mouse panel in Figure 4c, Extended Data Fig E8, and Extended Data
Fig E9, reproduced from released inputs only. Each clean script is a port or a
restyle of the original analysis script (the `R*.Rmd` files in this folder), kept
faithful to the published plotting style: manuscript palette, log10 axes with
10^n ticks, and the in-panel Wilcoxon / t-test brackets drawn by
`ggpubr::stat_compare_means`.

## How to run

```sh
# Public clone: reads released_data/ automatically.
Rscript reproduce/00_run_mouse.R

# Full private data: point at the complete folder, no code edits.
NUTRITION_DATA=/path/to/private/data Rscript reproduce/00_run_mouse.R

# Iterate on the CFU panels only, skipping the slow sequencing fits:
RUN_SEQ=false Rscript reproduce/00_run_mouse.R
```

Panels are written one PDF per panel into `results/`; expensive fits (MaAsLin2,
DESeq2, fgsea) cache into `intermediate_data/` and are reused on rerun.

### Packages

```r
renv::install(c("tidyverse","readxl","here","ggpubr","scales","ggrepel",
                "vegan","ashr","msigdbr"))
renv::install(c("bioc::Maaslin2","bioc::DESeq2","bioc::fgsea","bioc::edgeR",
                "bioc::org.Mm.eg.db","bioc::AnnotationDbi"))
renv::snapshot()
```

## Script organisation

```
reproduce/
  00_run_mouse.R            driver: helper, then 50-52, then (RUN_SEQ) 53-54
  50_figure_4c.R            F4c                       (style from R07 / R45)
  51_extdata_8.R            E8 g,h,j,k,l,m            (style from R07 / R08)
  52_extdata_9.R            E9 a,b,c,d,h              (style from R45 / R02)
  53_extdata_8_16s.R        E8 a,b,c,d,e,f  (16S)     (port of R01)
  54_extdata_9_rnaseq.R     E9 e,f,g        (RNA-seq) (port of R38)
  mouse/_mouse_helpers.R    data roots, palettes, log axis, theme, save, cache

  # Original publication scripts kept as reproduction reference:
  R01 16S (Fig S15) · R02 BMT 16S/CFU (S16F) · R07 no-fiber CFU (S16C)
  R08 chow (S16A-B) · R31 main CFU experiments · R38 GF RNA-seq · R45 smoothie
```

## Panel -> script -> output map

| Panel | Script | Released input | Output PDF | Reference |
|---|---|---|---|---|
| F4c upper | 50_figure_4c.R | workbook: Figure_4c_indiv_days_data | F4c_cfu_over_days.pdf | R31/R45 |
| F4c lower | 50_figure_4c.R | workbook: Figure_4c_trapezoidal_auc | F4c_auc.pdf | R31 |
| E8g | 51_extdata_8.R | workbook: Sup_Figure_8g_chow_consumed_day | E8g_chow_per_day.pdf | R08 |
| E8h | 51_extdata_8.R | workbook: Sup_Figure_8h_chow_consumed_auc | E8h_chow_auc.pdf | R08 |
| E8j | 51_extdata_8.R | workbook: Sup_Figure_8j_delayed_sucrose_a | E8j_delayed_sugar_auc.pdf | R31 |
| E8k | 51_extdata_8.R | workbook: Sup_Figure_8k_delayed_sucrose_C | E8k_cfu_over_days.pdf | R31 |
| E8l | 51_extdata_8.R | workbook: Sup_Figure_8l_no_fiber_chow_CFU | E8l_fiberfree_log10_cfu.pdf | R07 |
| E8m | 51_extdata_8.R | workbook: Sup_Figure_8m_mouse_percent_wei | E8m_weight_change.pdf | R31 |
| E8a | 53_extdata_8_16s.R | 204_mice_diet_healthy_all.csv | E8a_alpha_diversity.pdf | R01 |
| E8b | 53_extdata_8_16s.R | 204_mice_diet_healthy_all.csv | E8b_pcoa.pdf (+ E8b_permanova.txt) | R01 |
| E8c | 53_extdata_8_16s.R | 204_mice_diet_healthy_all.csv | E8c_abx_volcano_d3.pdf, _d6.pdf | R01 |
| E8e | 53_extdata_8_16s.R | 204_mice_diet_healthy_all.csv | E8e_sugar_volcano_d3.pdf, _d6.pdf | R01 |
| E8d | 53_extdata_8_16s.R | 204_mice_diet_healthy_all.csv | E8d_seqNNN.pdf (per ASV) | R01 |
| E8f | 53_extdata_8_16s.R | 204_mice_diet_healthy_all.csv | E8f_seqNNN.pdf (per ASV) | R01 |
| E9a | 52_extdata_9.R | workbook: Sup_Figure_9a_21_day_exp | E9a_21day_timecourse.pdf | R31 |
| E9b | 52_extdata_9.R | workbook: Sup_Figure_9b_smoothie | E9b_smoothie.pdf | R45 |
| E9c | 52_extdata_9.R | workbook: Sup_Figure_9a_21_day_exp (biapenem) | E9c_alternate_sugars.pdf | R31 |
| E9d time | 52_extdata_9.R | workbook: Sup_Figure_9d_time_course_monoc | E9d_monocolonization_timecourse.pdf | R31 |
| E9d AUC | 52_extdata_9.R | workbook: Sup_Figure_9d_auc_boxplot | E9d_monocolonization_auc.pdf | R31 |
| E9h | 52_extdata_9.R | workbook: Sup_Figure_9h_bmt_cfu (BM-only) | E9h_bmt_cfu.pdf | R02 |
| E9e | 54_extdata_9_rnaseq.R | raw_counts_matrix.tsv | E9e_gsea_lollipop.pdf | R38 |
| E9f | 54_extdata_9_rnaseq.R | raw_counts_matrix.tsv | E9f_pca.pdf | R38 |
| E9g | 54_extdata_9_rnaseq.R | raw_counts_matrix.tsv | E9g_deseq_volcano.pdf | R38 |

## Released data: purpose and key columns

**`Dai_mouse_figure_raw_data.xlsx`** — Nature source-data workbook, one sheet per
panel (sheet name = panel). Feeds all CFU, chow, weight and monocolonization
panels. Day is the trailing integer of the `Treatment + Day` label (parsed by
`add_day`). No PHI. (Full per-sheet column dictionary in the git history of this
file; the workbook sheets are self-describing via their headers.)

**`204_mice_diet_healthy_all.csv`** — long-format 16S table, the entry point for
all E8 sequencing panels (one ASV x sample per row).
- `Taxon` — `seqNNNN:taxonomy` ASV id with lineage; `seq` ids match the figure labels
- `sampleid` — stool sample id
- `relab` — relative abundance of the ASV in the sample
- `experiment_no`, `day` (0/3/6), `group`, `mouse_no`, `tube_no`, `investigator`
- `abx_treatment` (PBS/biapenem), `diet_treatment` (vehicle/sucrose)
- `simpson_reciprocal` — per-sample inverse Simpson alpha diversity (E8a)

**`raw_counts_matrix.tsv`** — RNA-seq gene-count matrix, entry point for E9 e/f/g.
- `Geneid` — versioned Ensembl mouse gene id (e.g. `ENSMUSG00000000001.5`)
- one column per sample, named `{exp}_{mouse}_{tissue}_IGO_...`; the sample sheet
  is derived from these names: gavage (mouse starting `1_`/`2_` = E. faecalis),
  diet (`1_`/`3_` = sucrose water), tissue (`_LI_` = large intestine = colon).
  All published panels use the LI samples.

## Analysis decisions carried from the originals

- 16S MaAsLin2: TSS + LOG, LM, `~ treatment * day + experiment_no`, min prevalence
  10%, min abundance 1e-4, q < 0.05. Two fits: antibiotic effect (PBS+vehicle
  reference) feeds E8c/E8d; sucrose effect (biapenem+vehicle reference) feeds
  E8e/E8f. Volcano significance: q < 0.05 and |coef| > 1.
- 16S PCoA: Bray-Curtis (`vegan::vegdist`) on the full ASV relab matrix, classical
  MDS; PERMANOVA (`adonis2`, 999 perms) on the biapenem arm per day, written to
  `results/E8b_permanova.txt`.
- RNA-seq DESeq2: design `~ treatment` on LI samples; `lfcShrink(type="ashr")`.
  Contrast 1 (E. faecalis vs PBS, regular water) feeds the E9e GSEA; contrast 2
  (sucrose vs regular, E. faecalis) feeds the E9g volcano and E9f PCA.
- GSEA (E9e): fgsea on MSigDB Hallmark (mouse), ranking metric
  `sign(log2FC) * -log10(pvalue)`, `set.seed(42)`, `nPermSimple = 10000`,
  size 15-500; lollipop of NES for pathways with padj < 0.05, dot size
  -log10(padj). 54 cross-checks its output against the shipped 9e/9g sheets.

## Notes for the run loop

- R is run by Angel in the sandbox; the session cannot execute R or the heavy
  Bioconductor packages, so treat 53/54 as needing a first run-and-report pass.
- 53/54 are faithful ports and keep the original magrittr `%>%` and spread/gather
  idioms; the CFU scripts (50-52) use the native pipe.
- E9a/E9c read one identical sheet (E9a = full antibiotic x sugar design, E9c =
  biapenem subset). The many-group time courses keep the palette and log axis but
  omit significance brackets; the originals report those tests in the text.
- E9h reproduces the BM-only (T-cell-depleted) recipients shown in the figure;
  the sheet also carries BM+Tcells, which is filtered out to match the panel.
- E8g released sheet dropped the per-day label, so E8g is the pooled boxplot
  rather than R08's per-day facet.
- Palettes: sucrose 4-group `gray76/#ffbcdc/gray32/deeppink2`; chow 2-group
  `gray32/deeppink2`; smoothie `gray76/bisque2/gray32/darkorange3`; volcano
  blue/black/red; lollipop royalblue/#E64B35.
```
