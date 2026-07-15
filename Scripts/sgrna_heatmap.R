# =============================================================================
# sgrna_heatmap.R
# -----------------------------------------------------------------------------
# Heatmap visualisation for CRISPR sgRNA screen results (MAGeCK output).
# Three heatmap types:
#   1. Manual candidate heatmap (ratio of two selected samples)
#   2. Full-sample gene-level heatmap (z-score, all samples)
#   3. Full-sample sgRNA-level heatmap (z-score, all samples)
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/sgrna-alignment-pipeline
# =============================================================================

source("sgrna_utils.R")

# =============================================================================
# PARAMETERS — edit here
# =============================================================================

DATA_DIR         <- "."
GENE_SUMMARY     <- file.path(DATA_DIR, "input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.gene_summary.txt")
SGRNA_SUMMARY    <- file.path(DATA_DIR, "input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.sgrna_summary.txt")
DATA_MATRIX      <- file.path(DATA_DIR, "input/sgRNAmatrix.tsv")
SAMPLES_FILE     <- file.path(DATA_DIR, "input/SampleAnnot.tsv")
PROJECT          <- "Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1"

REFERENCE        <- "D0_None_None"   # must have >= 1 expressed sgRNA
SCORE            <- 0.01             # FDR / p-value threshold
NB_MIN_SGRNA     <- 2                # min good sgRNAs per gene
SAMPLES_COMPARE  <- c("D8_Algox_FOLR1", "D8_UFO_FOLR1")  # for ratio heatmap

# =============================================================================
# SETUP
# =============================================================================

figs_base    <- file.path(DATA_DIR, "Figures/heatmaps", PROJECT)
results_base <- file.path(DATA_DIR, "results/heatmaps", PROJECT)
create_dirs(figs_base, results_base)

message("Loading data...")
counts   <- read_tsv(DATA_MATRIX)
annot    <- read_tsv(SAMPLES_FILE)
check_sample_names(counts, annot)

ordered_counts <- counts[order(rownames(counts)), rownames(annot)]
upm            <- normalise_upm(ordered_counts)

# DE gene filter (shared between heatmap 1 & 2)
genes    <- read_tsv(GENE_SUMMARY)
de_mask  <- (genes$neg.score < SCORE | genes$pos.score < SCORE) &
             (genes$neg.goodsgrna >= NB_MIN_SGRNA | genes$pos.goodsgrna >= NB_MIN_SGRNA)
de_genes <- rownames(genes)[de_mask]

# =============================================================================
# 1. MANUAL CANDIDATE HEATMAP (ratio of two samples)
# =============================================================================

message("[1/3] Manual candidate heatmap...")
gene_mat     <- upm[counts$Gene %in% de_genes, ]
gene_mat     <- gene_mat[gene_mat[, REFERENCE] > 0, ]

ratio_mat    <- (gene_mat[, SAMPLES_COMPARE[1]] + 0.001) /
                (gene_mat[, SAMPLES_COMPARE[2]] + 0.001)
ratio_mat    <- as.data.frame(ratio_mat)

col_ratio    <- colorRamp2(c(0, 1, 2), c("blue", "white", "red"))
ht_ratio     <- build_heatmap_sgrna(
  mat          = ratio_mat,
  row_labels   = counts[rownames(gene_mat), "Gene"],
  panel_name   = paste0(SAMPLES_COMPARE[1], "/", SAMPLES_COMPARE[2]),
  cluster_rows = FALSE,
  col_scale    = col_ratio
)
export_heatmap_pdf(ht_ratio,
  file.path(figs_base, paste0("Candidate_heatmap_Score_", SCORE, ".pdf")),
  n_rows = nrow(ratio_mat))
write_tsv(ratio_mat,
  file.path(results_base, paste0("Candidate_heatmap_Score_", SCORE, ".tsv")))

# =============================================================================
# 2. FULL-SAMPLE GENE HEATMAP (z-score)
# =============================================================================

message("[2/3] Full-sample gene heatmap...")
gene_mat2    <- upm[counts$Gene %in% de_genes, ]
gene_mat2    <- gene_mat2[gene_mat2[, REFERENCE] > 0, ]

ht_genes <- build_heatmap_sgrna(
  mat        = gene_mat2,
  row_labels = counts[rownames(gene_mat2), "Gene"],
  panel_name = "Genes_matrix"
)
export_heatmap_pdf(ht_genes,
  file.path(figs_base, paste0("Fullsample_Genes_heatmap_Score_", SCORE, ".pdf")),
  n_rows = nrow(gene_mat2))
write_tsv(as.data.frame(scale_rows(gene_mat2)),
  file.path(results_base, paste0("Fullsample_Genes_heatmap_Score_", SCORE, ".tsv")))

# =============================================================================
# 3. FULL-SAMPLE SGRNA HEATMAP (z-score)
# =============================================================================

message("[3/3] Full-sample sgRNA heatmap...")
sgrna        <- read_tsv(SGRNA_SUMMARY)
sgrna_mask   <- rownames(sgrna)[sgrna$p.twosided < SCORE]
sgrna_mat    <- upm[counts$Gene %in% sgrna_mask, ]
sgrna_mat    <- sgrna_mat[sgrna_mat[, REFERENCE] > 0, ]

ht_sgrna <- build_heatmap_sgrna(
  mat        = sgrna_mat,
  row_labels = counts[rownames(sgrna_mat), "Gene"],
  panel_name = "sgrna_matrix"
)
export_heatmap_pdf(ht_sgrna,
  file.path(figs_base, paste0("Fullsample_sgRNA_heatmap_Score_", SCORE, ".pdf")),
  n_rows = nrow(sgrna_mat))
write_tsv(as.data.frame(scale_rows(sgrna_mat)),
  file.path(results_base, paste0("Fullsample_sgRNA_heatmap_Score_", SCORE, ".tsv")))

message("Done. Output in: ", figs_base)
