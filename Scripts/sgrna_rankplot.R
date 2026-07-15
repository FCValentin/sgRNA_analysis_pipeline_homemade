# =============================================================================
# sgrna_rankplot.R
# -----------------------------------------------------------------------------
# Rank plot visualisation of MAGeCK gene-level enrichment scores.
# Produces negative and positive selection rank plots with top-gene labels.
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/sgrna-alignment-pipeline
# =============================================================================

source("sgrna_utils.R")

# =============================================================================
# PARAMETERS — edit here
# =============================================================================

DATA_DIR     <- "."
GENE_SUMMARY <- file.path(DATA_DIR, "input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.gene_summary.txt")
PROJECT      <- "Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1"
SCORE        <- 0.01    # FDR threshold for label display
TOP_RANK     <- 10      # only label genes with rank <= TOP_RANK

# =============================================================================
# SETUP
# =============================================================================

figs_dir <- file.path(DATA_DIR, "Figures/RankPlot", PROJECT)
create_dirs(figs_dir)

# =============================================================================
# RANK PLOTS
# =============================================================================

message("Loading gene summary...")
genes    <- read_tsv(GENE_SUMMARY)
pdf_path <- file.path(figs_dir, paste0("Genes_summary_Score_", SCORE, ".pdf"))
pdf(pdf_path, width = 10, height = 10)

# ── Negative selection rank plot ─────────────────────────────────────────────
neg_sig <- genes$neg.score < SCORE & genes$neg.rank <= TOP_RANK

plot(genes$neg.rank, -log10(genes$neg.score),
     type = "n", pch = 16, cex = 0.2,
     xlab = "Rank", ylab = "-log10(neg.score)",
     main = paste0("Negative selection rank — ", PROJECT))

if (any(neg_sig)) {
  text(genes$neg.rank[neg_sig],
       -log10(genes$neg.score[neg_sig]),
       labels = rownames(genes)[neg_sig],
       cex = 0.5, col = "green")
}
points(genes$neg.rank, -log10(genes$neg.score), cex = 0.2, col = "black", pch = 16)
abline(h = -log10(SCORE), lty = 2, col = "grey50")
abline(v = TOP_RANK,       lty = 2, col = "grey80")

# ── Positive selection rank plot ─────────────────────────────────────────────
pos_genes <- genes[order(genes$pos.rank), ]
pos_sig   <- pos_genes$pos.score < SCORE & pos_genes$pos.rank <= TOP_RANK

plot(pos_genes$pos.rank, -log10(pos_genes$pos.score),
     type = "n", pch = 16, cex = 0.2,
     xlab = "Rank", ylab = "-log10(pos.score)",
     main = paste0("Positive selection rank — ", PROJECT))

if (any(pos_sig)) {
  text(pos_genes$pos.rank[pos_sig],
       -log10(pos_genes$pos.score[pos_sig]),
       labels = rownames(pos_genes)[pos_sig],
       cex = 0.5, col = "red")
}
points(pos_genes$pos.rank, -log10(pos_genes$pos.score),
       cex = 0.2, col = "black", pch = 16)
abline(h = -log10(SCORE), lty = 2, col = "grey50")
abline(v = TOP_RANK,       lty = 2, col = "grey80")

dev.off()
message("Rank plot saved: ", pdf_path)
