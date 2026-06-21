# Mouse driver. Regenerates every mouse figure panel from released inputs:
# the source-data workbook (CFU, chow, weight, monocolonization), the long 16S
# relab table, and the RNA-seq count matrix.
#
# Run from the project root, with or without NUTRITION_DATA set:
#   Rscript reproduce/00_run_mouse.R
# Panels land in results/, one PDF per panel; cached fits in intermediate_data/.
#
# The CFU/behaviour scripts (50-52) are light and run anywhere. The sequencing
# scripts (53-54) need Bioconductor packages (Maaslin2, vegan; DESeq2, fgsea,
# org.Mm.eg.db, edgeR, ashr, msigdbr) and are the slow steps; set RUN_SEQ=FALSE
# in the environment to skip them while iterating on the CFU panels.

source(here::here("reproduce", "mouse", "_mouse_helpers.R"))

message("Reading mouse data from: ", nutrition_data_root())

source(here::here("reproduce", "50_figure_4c.R"))   # F4c
source(here::here("reproduce", "51_extdata_8.R"))   # E8 g,h,j,k,l,m
source(here::here("reproduce", "52_extdata_9.R"))   # E9 a,b,c,d,h

run_seq <- tolower(Sys.getenv("RUN_SEQ", unset = "true")) %in% c("true", "1", "yes")
if (run_seq) {
  source(here::here("reproduce", "53_extdata_8_16s.R"))     # E8 a,b,c,d,e,f (16S)
  source(here::here("reproduce", "54_extdata_9_rnaseq.R"))  # E9 e,f,g (RNA-seq)
} else {
  message("RUN_SEQ is off: skipping 16S (53) and RNA-seq (54) panels.")
}

message("Mouse panels complete. See results/.")
