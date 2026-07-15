# =============================================================================
# sgrna_volcano.R
# -----------------------------------------------------------------------------
# Volcano plot visualisation of MAGeCK gene-level and sgRNA-level results.
# Produces two plots per comparison:
#   1. Gene-level volcano (neg.lfc vs -log10(score))
#   2. sgRNA-level volcano (LFC vs -log10(FDR))
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/sgrna-alignment-pipeline
# =============================================================================

source("sgrna_utils.R")

# =============================================================================
# PARAMETERS — edit here
# =============================================================================

DATA_DIR      <- "."
GENE_SUMMARY  <- file.path(DATA_DIR, "input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.gene_summary.txt")
SGRNA_SUMMARY <- file.path(DATA_DIR, "input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.sgrna_summary.txt")
DATA_MATRIX   <- file.path(DATA_DIR, "input/sgRNAmatrix.tsv")
PROJECT       <- "Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1"
SCORE         <- 0.01    # FDR / p-value threshold
LFC_THRESHOLD <- 1.0     # absolute log2 FC threshold

# =============================================================================
# SETUP
# =============================================================================

figs_dir <- file.path(DATA_DIR, "Figures/Volcano", PROJECT)
create_dirs(figs_dir)

pdf_path <- file.path(figs_dir,
  paste0("VolcanoPlot_Score_", SCORE, "_FC_", LFC_THRESHOLD, ".pdf"))
pdf(pdf_path, width = 10, height = 10)

# =============================================================================
# 1. GENE-LEVEL VOLCANO
# =============================================================================

message("Plotting gene-level volcano...")
genes   <- read_tsv(GENE_SUMMARY)

pos_sig <- genes$pos.score < SCORE & genes$neg.lfc >  LFC_THRESHOLD
neg_sig <- genes$neg.score < SCORE & genes$neg.lfc < -LFC_THRESHOLD

plot(genes$neg.lfc, -log10(genes$neg.score),
     type = "n",
     xlab = "neg.lfc",
     ylab = "-log10(neg.score)",
     main = paste0("Gene-level volcano — ", PROJECT))
abline(h = -log10(SCORE),    lty = 2, col = "grey50")
abline(v = c(-LFC_THRESHOLD, LFC_THRESHOLD), lty = 2, col = "grey70")

if (any(pos_sig))
  text(genes$neg.lfc[pos_sig], -log10(genes$pos.score[pos_sig]),
       labels = rownames(genes)[pos_sig], cex = 0.2, col = "red")
if (any(neg_sig))
  text(genes$neg.lfc[neg_sig], -log10(genes$neg.score[neg_sig]),
       labels = rownames(genes)[neg_sig], cex = 0.5, col = "green")

points(c(genes$neg.lfc, genes$pos.lfc),
       c(-log10(genes$neg.score), -log10(genes$pos.score)),
       cex = 0.2, col = "black", pch = 16)

# =============================================================================
# 2. SGRNA-LEVEL VOLCANO
# =============================================================================

message("Plotting sgRNA-level volcano...")
sgrna  <- read_tsv(SGRNA_SUMMARY)
counts <- read_tsv(DATA_MATRIX)

sgrna_pos <- sgrna$FDR < SCORE & sgrna$LFC >  LFC_THRESHOLD
sgrna_neg <- sgrna$FDR < SCORE & sgrna$LFC < -LFC_THRESHOLD

plot(sgrna$LFC, -log10(sgrna$FDR),
     type = "n",
     xlab = "LFC",
     ylab = "-log10(FDR)",
     main = paste0("sgRNA-level volcano — ", PROJECT))
abline(h = -log10(SCORE),    lty = 2, col = "grey50")
abline(v = c(-LFC_THRESHOLD, LFC_THRESHOLD), lty = 2, col = "grey70")

if (any(sgrna_pos))
  text(sgrna$LFC[sgrna_pos], -log10(sgrna$FDR[sgrna_pos]),
       labels = counts[rownames(sgrna)[sgrna_pos], "Gene"],
       cex = 0.3, col = "red")
if (any(sgrna_neg))
  text(sgrna$LFC[sgrna_neg], -log10(sgrna$FDR[sgrna_neg]),
       labels = counts[rownames(sgrna)[sgrna_neg], "Gene"],
       cex = 0.3, col = "green")

points(sgrna$LFC, -log10(sgrna$FDR), cex = 0.2, col = "black", pch = 16)

dev.off()
message("Volcano plots saved: ", pdf_path)
