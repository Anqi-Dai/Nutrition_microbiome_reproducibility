# Extended Data Figure E8, 16S panels.
#
# Faithful port of R01 (Figure S15 in the original numbering). The long-format
# relab table was split for the workbook into two consolidated sheets, which this
# script joins back into the single `healthyall` table R01 worked from:
#   mouse_16s_asv_relab    wide ASV relative abundance (Taxon rows x sample cols)
#   mouse_16s_sample_meta  one row per sample (treatment, day, alpha diversity, ...)
#
# One table feeds six panels:
#   E8a  alpha diversity (inverse Simpson), paired day 1 vs day 3 Wilcoxon
#   E8b  PCoA of Bray-Curtis distances, faceted by day (+ PERMANOVA to a log)
#   E8c  MaAsLin2 volcano, antibiotic effect (PBS+vehicle reference), d3 & d6
#   E8e  MaAsLin2 volcano, sucrose-on-top-of-antibiotic effect, d3 & d6
#   E8d  relative-abundance trajectories of the antibiotic-enriched ASVs
#   E8f  relative-abundance trajectories of the sucrose-enriched ASVs
#
# Days read 1/3/6 (the original day-0 baseline label was a wet-lab error, now
# corrected in mouse_16s_sample_meta). alpha, PCoA and relab are deterministic
# transforms; only the two MaAsLin2 fits are expensive, so they cache under
# intermediate_data and are reused on rerun.
#
# Ported from the original, so this keeps the magrittr pipe and the spread/gather
# idioms it relied on. Heavy dependencies (Maaslin2, vegan) install via:
#   renv::install(c("vegan","ggrepel")); renv::install("bioc::Maaslin2")

source(here::here("reproduce", "mouse", "_mouse_helpers.R"))
suppressPackageStartupMessages({
  library(tibble)
  library(purrr)      # map / iwalk over the per-day and per-ASV splits
  library(vegan)
  library(ggrepel)
  library(Maaslin2)
})

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)

# Reconstruct R01's long table from the two consolidated 16S sheets. Days are
# 1/3/6 in the source data (the original day-0 baseline label was a wet-lab error,
# corrected to day 1 in mouse_16s_sample_meta), so no relabel is needed here.
relab_wide <- read_mouse_sheet("mouse_16s_asv_relab")
meta_sheet <- read_mouse_sheet("mouse_16s_sample_meta")
healthyall <- relab_wide %>%
  pivot_longer(-Taxon, names_to = "sampleid", values_to = "relab") %>%
  left_join(meta_sheet, by = "sampleid")

meta <- healthyall %>%
  select(-Taxon, -relab) %>%
  distinct()

grp_levels <- c("PBS__vehicle", "PBS__sucrose", "biapenem__vehicle", "biapenem__sucrose")

# E8a: alpha diversity --------------------------------------------------------
alpha_df <- meta %>%
  mutate(day = factor(day), experiment_no = factor(experiment_no)) %>%
  arrange(abx_treatment, desc(diet_treatment), day) %>%
  mutate(xvar = str_glue("{abx_treatment}__{diet_treatment}__{day}"),
         grp  = str_glue("{abx_treatment}__{diet_treatment}"),
         grp  = factor(grp, levels = grp_levels))

alpha_change <- alpha_df %>%
  ggboxplot(x = "xvar", y = "simpson_reciprocal", add = "jitter",
            xlab = "", ylab = "alpha diversity (Simpson reciprocal)",
            add.params = list(shape = "experiment_no", alpha = 0.8),
            width = 0.6, color = "grp") +
  scale_color_manual(values = pal_sucrose4) +
  stat_compare_means(comparisons = list(
    c("biapenem__vehicle__1", "biapenem__vehicle__3"),
    c("biapenem__sucrose__1", "biapenem__sucrose__3"),
    c("PBS__vehicle__1", "PBS__vehicle__3"),
    c("PBS__sucrose__1", "PBS__sucrose__3")),
    paired = TRUE, label = "p.signif", method = "wilcox.test",
    tip.length = 0.04, exact = TRUE, correct = TRUE) +
  scale_x_discrete(labels = rep(c(1, 3, 6), 4)) +
  scale_y_sqrt() +
  theme_light() +
  theme(legend.position = "none", aspect.ratio = 1 / 1.3)

save_panel(alpha_change, "E8a_alpha_diversity.pdf", width = 5.0, height = 3.4)

# E8b: PCoA on Bray-Curtis ----------------------------------------------------
relab_mat <- healthyall %>%
  select(Taxon, sampleid, relab) %>%
  spread(key = "Taxon", value = "relab", fill = 0) %>%
  column_to_rownames("sampleid")

dist_bc <- vegdist(relab_mat, method = "bray")
eigen <- cmdscale(dist_bc, eig = TRUE)$eig
percent_var <- signif(eigen / sum(eigen), 3) * 100
bc <- cmdscale(dist_bc, k = 2)

beta_all <- bc %>%
  as.data.frame() %>%
  rownames_to_column("sampleid") %>%
  full_join(meta, by = "sampleid") %>%
  mutate(day = factor(day), experiment_no = factor(experiment_no)) %>%
  arrange(abx_treatment, desc(diet_treatment), day) %>%
  mutate(grp = str_glue("{abx_treatment}__{diet_treatment}"),
         grp = factor(grp, levels = grp_levels),
         experiment_no = str_glue("Expe. {experiment_no}"),
         day = str_glue("D{day}"))

pcoa <- beta_all %>%
  ggscatter(x = "V1", y = "V2", color = "grp", shape = "experiment_no", alpha = 0.8) +
  scale_color_manual(values = pal_sucrose4) +
  facet_grid(. ~ day) +
  xlab(paste0("PC 1 [", percent_var[1], "%]")) +
  ylab(paste0("PC 2 [", percent_var[2], "%]")) +
  theme_light() +
  theme(aspect.ratio = 1, legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

save_panel(pcoa, "E8b_pcoa.pdf", width = 6.5, height = 3.0)

# PERMANOVA (biapenem arm, per day), written to a log for the record.
permanova <- healthyall %>%
  select(Taxon, sampleid, relab, abx_treatment, diet_treatment, experiment_no, day) %>%
  spread(key = "Taxon", value = "relab", fill = 0) %>%
  filter(abx_treatment == "biapenem") %>%
  column_to_rownames("sampleid") %>%
  split(.$day) %>%
  map(function(df) {
    d <- vegdist(df %>% select(-c(abx_treatment, diet_treatment, experiment_no, day)),
                 method = "bray")
    meta_d <- df %>% select(diet_treatment, experiment_no)
    adonis2(d ~ diet_treatment + experiment_no, data = meta_d, permutations = 999)
  })
capture.output(print(permanova), file = cache_path("E8b_permanova.txt"))

# MaAsLin2 helper: fit once, cache the output dir, return all_results.tsv -------
run_maaslin <- function(asv_long, treatment_levels, out_name, fixed_effects) {
  out_dir <- cache_path(out_name)
  all_res <- file.path(out_dir, "all_results.tsv")
  if (!file.exists(all_res)) {
    df_t <- asv_long %>%
      mutate(treatment = factor(str_glue("{abx_treatment}_{diet_treatment}"),
                                levels = treatment_levels))
    n_samp <- df_t %>% distinct(sampleid) %>% nrow()
    kept <- df_t %>%
      group_by(Taxon) %>%
      summarise(perc = round(sum(relab > 0.0001) / n_samp * 100, 0), .groups = "drop") %>%
      filter(perc > 10) %>% pull(Taxon)
    wide <- df_t %>%
      filter(Taxon %in% kept) %>%
      select(Taxon, sampleid, relab, treatment, day, experiment_no) %>%
      spread("Taxon", "relab", fill = 0) %>%
      mutate(day = factor(day), experiment_no = factor(experiment_no))

    md <- wide %>% select(sampleid, treatment, day, experiment_no) %>% column_to_rownames("sampleid")
    design <- model.matrix(~ treatment * day + experiment_no, data = md)
    md <- cbind(md, design)
    feats <- wide %>% select(sampleid, starts_with("seq")) %>% column_to_rownames("sampleid")

    Maaslin2(input_data = feats, input_metadata = md, output = out_dir,
             fixed_effects = fixed_effects, random_effects = NULL,
             normalization = "TSS", transform = "LOG", analysis_method = "LM",
             max_significance = 0.05, min_abundance = 0.0001, min_prevalence = 0.10,
             plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1)
  }
  read_tsv(all_res, show_col_types = FALSE) %>%
    separate(feature, into = c("seq_id", "taxa"), sep = "\\.", remove = FALSE, extra = "merge")
}

asv_long <- healthyall %>% select(Taxon, sampleid, relab, abx_treatment, diet_treatment, day, experiment_no)

# E8c: antibiotic effect (PBS+vehicle reference) ------------------------------
abx_fixed <- c("treatmentPBS_sucrose", "treatmentbiapenem_vehicle", "treatmentbiapenem_sucrose",
               "day3", "day6", "experiment_no2",
               "treatmentPBS_sucrose:day3", "treatmentbiapenem_vehicle:day3", "treatmentbiapenem_sucrose:day3",
               "treatmentPBS_sucrose:day6", "treatmentbiapenem_vehicle:day6", "treatmentbiapenem_sucrose:day6")
all_res_abx <- run_maaslin(asv_long,
                           c("PBS_vehicle", "PBS_sucrose", "biapenem_vehicle", "biapenem_sucrose"),
                           "R01_maaslin2_abx_effect", abx_fixed)

# E8e: sucrose-on-top-of-antibiotic effect (biapenem+vehicle reference) --------
sugar_fixed <- c("treatmentbiapenem_sucrose", "treatmentPBS_vehicle", "treatmentPBS_sucrose",
                 "day3", "day6", "experiment_no2",
                 "treatmentbiapenem_sucrose:day3", "treatmentPBS_vehicle:day3", "treatmentPBS_sucrose:day3",
                 "treatmentbiapenem_sucrose:day6", "treatmentPBS_vehicle:day6", "treatmentPBS_sucrose:day6")
all_res_sugar <- run_maaslin(asv_long,
                             c("biapenem_vehicle", "biapenem_sucrose", "PBS_vehicle", "PBS_sucrose"),
                             "R01_maaslin2_sugar_effect", sugar_fixed)

# Volcano builder, shared by E8c and E8e --------------------------------------
volcano <- function(all_res, contrast, title) {
  all_res %>%
    filter(metadata == contrast) %>%
    mutate(sig_ = if_else(qval < 0.05 & coef > 1, "pos",
                  if_else(qval < 0.05 & coef < -1, "neg", "not")),
           neglog10_q = -log10(qval)) %>%
    ggscatter(x = "coef", y = "neglog10_q", color = "sig_", alpha = 0.5, shape = 16,
              xlab = "MaAsLin2 model coefficient", ylab = "- log10(q value)", title = title) +
    scale_color_manual(values = pal_volcano) +
    geom_vline(xintercept = 0, color = "black") +
    geom_vline(xintercept = c(1, -1), color = "gray", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "gray", linetype = "dashed") +
    geom_text_repel(aes(label = ifelse(qval < 0.05 & coef > 1, seq_id, "")), color = "red", seed = 1) +
    geom_text_repel(aes(label = ifelse(qval < 0.05 & coef < -1, seq_id, "")), color = "blue", seed = 1) +
    theme(aspect.ratio = 1 / 2, legend.position = "none")
}

save_panel(volcano(all_res_abx, "treatmentbiapenem_vehicle:day3", "Day 1 vs 3"),
           "E8c_abx_volcano_d3.pdf", width = 4.2, height = 3.0)
save_panel(volcano(all_res_abx, "treatmentbiapenem_vehicle:day6", "Day 1 vs 6"),
           "E8c_abx_volcano_d6.pdf", width = 4.2, height = 3.0)
save_panel(volcano(all_res_sugar, "treatmentbiapenem_sucrose:day3", "Day 1 vs 3"),
           "E8e_sugar_volcano_d3.pdf", width = 4.2, height = 3.0)
save_panel(volcano(all_res_sugar, "treatmentbiapenem_sucrose:day6", "Day 1 vs 6"),
           "E8e_sugar_volcano_d6.pdf", width = 4.2, height = 3.0)

# E8d / E8f: relative-abundance trajectories for the highlighted ASVs ----------
relab_trajectories <- function(all_res, contrasts, hold_var, hold_levels, prefix) {
  feats <- all_res %>%
    filter(metadata %in% contrasts, qval < 0.05, abs(coef) > 1) %>%
    distinct(feature) %>% pull(feature)
  if (length(feats) == 0) return(invisible(NULL))

  facet_col <- if (hold_var == "diet_treatment") "abx_treatment" else "diet_treatment"
  dat <- healthyall %>%
    mutate(feature = str_replace_all(str_replace(Taxon, "\\:", "."), "\\;", ".")) %>%
    filter(feature %in% feats)

  dat %>%
    split(.$feature) %>%
    iwalk(function(df, feat) {
      seqid <- str_extract(feat, "^seq\\d+")
      p <- df %>%
        filter(.data[[hold_var]] %in% hold_levels) %>%
        mutate(facet_var = factor(.data[[facet_col]]),
               u_mouse_id = str_glue("{experiment_no}_{mouse_no}")) %>%
        ggboxplot(x = "day", y = "relab", outlier.shape = NA, xlab = "Day",
                  ylab = "Relative abundance", title = seqid, width = 0.6) +
        geom_point(aes(group = u_mouse_id), alpha = 0.3, shape = 16) +
        geom_line(aes(group = u_mouse_id), alpha = 0.3) +
        facet_grid(. ~ facet_var) +
        theme(legend.position = "none", aspect.ratio = 1.5)
      save_panel(p, str_glue("{prefix}_{seqid}.pdf"), width = 3.2, height = 3.0)
    })
}

# E8d: antibiotic-enriched taxa, held in the vehicle arm, faceted by antibiotic.
relab_trajectories(all_res_abx,
                   c("treatmentbiapenem_vehicle:day3", "treatmentbiapenem_vehicle:day6"),
                   hold_var = "diet_treatment", hold_levels = "vehicle", prefix = "E8d")

# E8f: sucrose-enriched taxa, held in the biapenem arm, faceted by diet.
relab_trajectories(all_res_sugar,
                   c("treatmentbiapenem_sucrose:day3", "treatmentbiapenem_sucrose:day6"),
                   hold_var = "abx_treatment", hold_levels = "biapenem", prefix = "E8f")

highlighted_seqs <- function(all_res, contrasts) {
  all_res %>%
    filter(metadata %in% contrasts, qval < 0.05, abs(coef) > 1) %>%
    distinct(seq_id) %>% pull(seq_id) %>% sort()
}
message("E8d antibiotic-enriched ASVs: ",
        paste(highlighted_seqs(all_res_abx,
              c("treatmentbiapenem_vehicle:day3", "treatmentbiapenem_vehicle:day6")), collapse = ", "))
message("E8f sucrose-enriched ASVs: ",
        paste(highlighted_seqs(all_res_sugar,
              c("treatmentbiapenem_sucrose:day3", "treatmentbiapenem_sucrose:day6")), collapse = ", "))
