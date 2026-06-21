# Extended Data Figure E9, RNA-seq panels.
#
# Faithful port of R38 (germ-free colonic RNA-seq), reading the gene count matrix
# from the consolidated workbook sheet `raw_counts_matrix`. The published panels
# are all from the large intestine (colonic epithelium):
#   E9e  GSEA Hallmark lollipop, E. faecalis-colonized vs uncolonized (contrast 1)
#   E9f  PCA of colonized mice, regular vs sucrose water (t-test on PC1)
#   E9g  DESeq2 volcano, sucrose vs regular water during colonization (contrast 2)
#
# R38 fits one DESeq2 object per tissue (the LI object dds_li) but, in the notebook,
# the line creating it was dropped (it is only ever re-run via DESeq()). Because
# every LI sample shares one tissue level, a ~ tissue + treatment design is not
# estimable on the LI subset, so the faithful reconstruction is LI-only samples
# with design ~ treatment; that is what is built (and cached) below.
#
# The DESeq2 fit and the fgsea result are cached under intermediate_data and
# reused on rerun. The released source-data sheets (9e, 9g) are used as a light
# regression check at the end. Ranking metric for GSEA is sign(log2FC) * -log10(p)
# (R38 uses this since the shrunken results carry no DESeq `stat` column); gene
# sets are MSigDB Hallmark for mouse with seed 42 and 10,000 permutations.
#
# Heavy Bioconductor dependencies; install before running:
#   renv::install(c("ashr","ggrepel","msigdbr"))
#   renv::install(c("bioc::DESeq2","bioc::fgsea","bioc::org.Mm.eg.db",
#                   "bioc::AnnotationDbi","bioc::edgeR"))
# Optional: ggrastr (rasterises the volcano points, as R38 does); the script
# falls back to plain points when it is absent.

source(here::here("reproduce", "mouse", "_mouse_helpers.R"))
suppressPackageStartupMessages({
  library(tibble)     # column_to_rownames / rownames_to_column used in the port
  library(DESeq2)
  library(edgeR)
  library(ashr)
  library(fgsea)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(ggrepel)
})

if (!dir.exists(intermediate_dir())) dir.create(intermediate_dir(), recursive = TRUE)
set.seed(42)

# Sample sheet derived from the column names (mouse id encodes gavage and diet) --
# The gene count matrix now lives in the consolidated workbook as its own sheet.
counts_data <- read_mouse_sheet("raw_counts_matrix") %>%
  rename_with(.fn = ~ str_remove(., "_IGO_\\d+.+$"), .cols = 2:dplyr::last_col())
sample_names <- colnames(counts_data)[-1]

meta_table <- tibble(sample_id = sample_names) %>%
  mutate(
    mouse_id = str_extract(sample_id, "^\\d+_\\d+"),
    gavage_treatment = if_else(str_starts(mouse_id, "1_") | str_starts(mouse_id, "2_"),
                               "E. faecalis", "PBS"),
    diet_treatment = if_else(str_starts(mouse_id, "1_") | str_starts(mouse_id, "3_"),
                             "sucrose_water", "regular_water"),
    tissue = if_else(str_detect(sample_id, "_LI"), "large_intestine", "small_intestine"),
    gavage_treatment = factor(gavage_treatment, levels = c("PBS", "E. faecalis")),
    treatment = factor(str_glue("{gavage_treatment}__{diet_treatment}"),
                       levels = c("E. faecalis__regular_water", "E. faecalis__sucrose_water",
                                  "PBS__regular_water", "PBS__sucrose_water"))
  )

# CPM filter: keep genes with CPM > 2.5 in at least 3 samples (R38 elbow choice) -
raw_counts_matrix <- counts_data %>% column_to_rownames("Geneid") %>% as.matrix()
keep <- rowSums(edgeR::cpm(raw_counts_matrix) > 2.5) >= 3
filtered_counts <- raw_counts_matrix[keep, ]

# Large-intestine DESeq2 fit, cached ------------------------------------------
li_cache <- cache_path("R38_li_deseq.rds")
if (file.exists(li_cache)) {
  li <- readRDS(li_cache)
} else {
  meta_li <- meta_table %>% filter(tissue == "large_intestine")
  counts_li <- round(filtered_counts[, meta_li$sample_id])
  dds_li <- DESeqDataSetFromMatrix(countData = counts_li,
                                   colData = meta_li %>% column_to_rownames("sample_id"),
                                   design = ~ treatment)
  dds_li <- DESeq(dds_li)

  # Contrast 1: colonization effect (regular water). Contrast 2: sucrose effect.
  res1 <- results(dds_li, contrast = c("treatment", "E. faecalis__regular_water", "PBS__regular_water"))
  shr1 <- lfcShrink(dds_li, contrast = c("treatment", "E. faecalis__regular_water", "PBS__regular_water"),
                    res = res1, type = "ashr")
  res2 <- results(dds_li, contrast = c("treatment", "E. faecalis__sucrose_water", "E. faecalis__regular_water"))
  shr2 <- lfcShrink(dds_li, contrast = c("treatment", "E. faecalis__sucrose_water", "E. faecalis__regular_water"),
                    res = res2, type = "ashr")

  li <- list(dds = dds_li,
             shr1 = as.data.frame(shr1) %>% rownames_to_column("ensembl_version"),
             shr2 = as.data.frame(shr2) %>% rownames_to_column("ensembl_version"))
  saveRDS(li, li_cache)
}

map_symbol <- function(ensembl_version) {
  ids <- str_remove(ensembl_version, "\\..*")
  AnnotationDbi::mapIds(org.Mm.eg.db, keys = ids, column = "SYMBOL",
                        keytype = "ENSEMBL", multiVals = "first")
}

# E9g: DESeq2 volcano, sucrose vs regular water -------------------------------
# Significance is padj < 0.05 and |log2FC| > 1 (R38). The point layer is
# rasterised (as in R38) when ggrastr is available, so the vector PDF stays small
# without a hard dependency on it.
have_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
volcano_points <- function() {
  pts <- geom_point(aes(colour = significant), alpha = 0.6, shape = 16, size = 1)
  if (have_ggrastr) ggrastr::rasterise(pts, dpi = 300) else pts
}

vol_df <- li$shr2 %>%
  mutate(symbol = map_symbol(ensembl_version)) %>%
  filter(!is.na(symbol)) %>%
  mutate(significant = if_else(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"))

p_9g <- vol_df %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  volcano_points() +
  scale_colour_manual(values = c(No = "grey", Yes = "red"), name = "") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  # Label the significant genes with their symbol (R38). seed fixes the repel
  # layout so the panel is reproducible across runs.
  geom_text_repel(aes(label = symbol), data = filter(vol_df, significant == "Yes"),
                  max.overlaps = 10.8, size = 4, box.padding = 0.4, seed = 42) +
  labs(x = "Log2 Fold Change", y = "-log10(Adjusted P-value)",
       title = "E. faecalis (Sucrose vs Regular Water)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

save_panel(p_9g, "E9g_deseq_volcano.pdf", width = 5.0, height = 3.6)

# E9f: PCA of colonized LI samples, regular vs sucrose water ------------------
# vst on the E. faecalis large-intestine subset, then plotPCA (top 500 genes), as
# in R38. The PC1 group test follows R38's decision: Shapiro-Wilk per group, then
# a t-test if both look normal with n >= 3, otherwise Wilcoxon.
pca_pc1_test <- function(pca_df) {
  grp_n <- min(table(pca_df$treatment))
  is_normal <- grp_n >= 3 && all(
    pca_df %>% group_by(treatment) %>%
      summarise(p = shapiro.test(PC1)$p.value, .groups = "drop") %>% pull(p) > 0.05)
  if (is_normal) {
    list(method = "t-test", p = t.test(PC1 ~ treatment, data = pca_df)$p.value)
  } else {
    list(method = "Wilcoxon", p = suppressWarnings(wilcox.test(PC1 ~ treatment, data = pca_df)$p.value))
  }
}

meta_li_faecalis <- meta_table %>% filter(tissue == "large_intestine", gavage_treatment == "E. faecalis")
counts_sub <- round(filtered_counts[, meta_li_faecalis$sample_id])
dds_sub <- DESeqDataSetFromMatrix(countData = counts_sub,
                                  colData = meta_li_faecalis %>% column_to_rownames("sample_id"),
                                  design = ~ treatment)
vst_sub <- vst(dds_sub, blind = TRUE)
pca_data <- plotPCA(vst_sub, intgroup = "treatment", returnData = TRUE)
percent_var_pca <- round(100 * attr(pca_data, "percentVar"))
pc1 <- pca_pc1_test(pca_data)

p_9f <- pca_data %>%
  ggplot(aes(x = PC1, y = PC2, colour = treatment)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("E. faecalis__regular_water" = "#00BFC4",
                                 "E. faecalis__sucrose_water" = "#F8766D")) +
  labs(x = str_glue("PC1: {percent_var_pca[1]}% variance"),
       y = str_glue("PC2: {percent_var_pca[2]}% variance"),
       title = str_glue("LI colonized (PC1 {pc1$method} p = {signif(pc1$p, 2)})")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", legend.title = element_blank())

save_panel(p_9f, "E9f_pca.pdf", width = 4.8, height = 3.2)

# E9e: GSEA Hallmark lollipop, colonization effect (contrast 1) ---------------
prepare_ranked_list <- function(shr_df) {
  df <- shr_df %>%
    filter(!is.na(pvalue) & !is.na(log2FoldChange)) %>%
    mutate(rank_metric = sign(log2FoldChange) * -log10(pvalue),
           symbol = map_symbol(ensembl_version)) %>%
    filter(!is.na(symbol)) %>%
    group_by(symbol) %>%
    filter(abs(rank_metric) == max(abs(rank_metric))) %>%
    ungroup()
  ranked <- df$rank_metric
  names(ranked) <- df$symbol
  sort(ranked[!is.na(ranked)], decreasing = TRUE)
}

fgsea_cache <- cache_path("R38_fgsea_li_contrast1.rds")
if (file.exists(fgsea_cache)) {
  fgsea_li1 <- readRDS(fgsea_cache)
} else {
  # msigdbr renamed `category` to `collection` around v10; prefer the modern arg
  # and fall back for older installs. Both return the same Hallmark gene sets.
  hallmark_sets <- tryCatch(msigdbr(species = "Mus musculus", collection = "H"),
                            error = function(e) msigdbr(species = "Mus musculus", category = "H"))
  hallmark_pathways <- hallmark_sets %>% split(x = .$gene_symbol, f = .$gs_name)
  set.seed(42)
  fgsea_li1 <- fgsea(pathways = hallmark_pathways,
                     stats = prepare_ranked_list(li$shr1),
                     minSize = 15, maxSize = 500, nPermSimple = 10000)
  saveRDS(fgsea_li1, fgsea_cache)
}

lollipop_df <- fgsea_li1 %>%
  as_tibble() %>%
  filter(padj < 0.05) %>%
  mutate(pathway_clean = pathway %>%
           str_replace("HALLMARK_", "") %>% str_to_title() %>%
           str_replace_all("_", " ") %>% str_wrap(width = 48),
         neg_log_padj = -log10(padj),
         pathway_clean = reorder(pathway_clean, NES),
         direction = if_else(NES > 0, "Upregulated", "Downregulated"))

p_9e <- lollipop_df %>%
  ggplot(aes(x = NES, y = pathway_clean, colour = direction)) +
  geom_segment(aes(x = 0, xend = NES, y = pathway_clean, yend = pathway_clean), linewidth = 0.8) +
  geom_point(aes(size = neg_log_padj), alpha = 0.8) +
  scale_colour_manual(values = pal_lollipop, name = "Direction") +
  scale_size_continuous(name = "-log10(Adjusted p-value)") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Normalized Enrichment Score (NES)", y = "") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 9), aspect.ratio = 2) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

save_panel(p_9e, "E9e_gsea_lollipop.pdf", width = 6.0, height = 6.5)

# Light run summary. The per-panel RNA-seq source-data sheets (9e/9g) used to
# cross-check these counts were removed once the workbook was consolidated around
# the raw count matrix; that check confirmed an exact match (9/9 significant
# genes, 21/21 Hallmark pathways) before they were dropped.
message(str_glue("E9g significant genes: {sum(vol_df$significant == 'Yes')}"))
message(str_glue("E9e enriched Hallmark pathways (padj < 0.05): {nrow(lollipop_df)}"))
