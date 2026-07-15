# =============================================================================
# sgrna_mle.R
# -----------------------------------------------------------------------------
# MLE (Maximum Likelihood Estimation) analysis of MAGeCK MLE output.
# For each design matrix variable:
#   - Beta-score filtered gene heatmap
#   - Beta-score distribution histogram
#   - GO enrichment for up/down/all genes
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
DATA_MATRIX   <- file.path(DATA_DIR, "input/sgRNAmatrix.tsv")
GENE_SUMMARY  <- file.path(DATA_DIR, "input/test.mle.gene_summary.txt")
SGRNA_SUMMARY <- file.path(DATA_DIR, "input/test.mle.sgrna_summary.txt")
SAMPLES_FILE  <- file.path(DATA_DIR, "input/SampleAnnot_binary.tsv")
ORTHOLOG_FILE <- file.path(DATA_DIR, "input/HumanOrthologs.tsv")
PROJECT       <- "MLE"

REFERENCE     <- "D0_None_None"   # reference sample
BETA_SCORE    <- 2                # |beta| threshold

# =============================================================================
# SETUP
# =============================================================================

figs_dir    <- file.path(DATA_DIR, "Figures/MLE")
results_dir <- file.path(DATA_DIR, "results/MLE")
create_dirs(figs_dir, results_dir)

message("Loading data...")
counts    <- read_tsv(DATA_MATRIX)
annot     <- read_tsv(SAMPLES_FILE)
genes     <- read_tsv(GENE_SUMMARY)
sgrna     <- read.table(SGRNA_SUMMARY, row.names = 2, header = TRUE)
orthologs <- read_tsv(ORTHOLOG_FILE)

check_sample_names(counts, annot)
ordered_counts <- counts[order(rownames(counts)), rownames(annot)]
upm            <- normalise_upm(ordered_counts)

# =============================================================================
# LOOP OVER DESIGN MATRIX VARIABLES (skip intercept col 1)
# =============================================================================

for (i in seq(2, ncol(annot))) {
  var_name <- colnames(annot)[i]
  beta_col <- paste0(var_name, ".beta")
  message("\n[", i - 1, "/", ncol(annot) - 1, "] Variable: ", var_name)

  if (!beta_col %in% colnames(genes)) {
    warning("Beta column not found: ", beta_col, " — skipping.")
    next
  }

  # Filter genes by beta score
  de_mask  <- abs(genes[, beta_col]) > BETA_SCORE
  de_genes <- rownames(genes)[de_mask]
  gene_mat <- upm[counts$Gene %in% de_genes, ]
  gene_mat <- gene_mat[gene_mat[, REFERENCE] > 0, ]

  if (nrow(gene_mat) == 0) {
    message("  No genes pass beta score threshold — skipping heatmap.")
    next
  }

  # Heatmap + histogram
  ht <- build_heatmap_sgrna(
    mat        = gene_mat,
    row_labels = counts[rownames(gene_mat), "Gene"],
    panel_name = "Genes_matrix"
  )
  pdf_path <- file.path(figs_dir,
    paste0(var_name, "_Genes_heatmap_BetaScore_", BETA_SCORE, ".pdf"))
  pdf(pdf_path, width = 10, height = pdf_height_from_nrow(nrow(gene_mat)))
  hist(genes[, beta_col], xlab = "Beta score",
       main = paste("Beta score distribution —", var_name))
  print(ht)
  dev.off()
  message("  Heatmap saved: ", pdf_path)

  write_tsv(as.data.frame(scale_rows(gene_mat)),
    file.path(results_dir,
      paste0(var_name, "_Genes_heatmap_BetaScore_", BETA_SCORE, ".tsv")))

  # GO enrichment
  up_genes   <- rownames(genes)[!is.na(genes[, beta_col]) & genes[, beta_col] >  BETA_SCORE]
  down_genes <- rownames(genes)[!is.na(genes[, beta_col]) & genes[, beta_col] < -BETA_SCORE]
  up_genes   <- up_genes[up_genes %in% orthologs$Symbol]
  down_genes <- down_genes[down_genes %in% orthologs$Symbol]
  up_ortho   <- orthologs[orthologs$Symbol %in% up_genes, ]
  down_ortho <- orthologs[orthologs$Symbol %in% down_genes, ]

  write_tsv(up_ortho,   file.path(results_dir, paste0(var_name, "_List_UpOrthologs.tsv")))
  write_tsv(down_ortho, file.path(results_dir, paste0(var_name, "_List_DownOrthologs.tsv")))

  run_gprofiler(up_ortho$Gene.Group.Identifier,
    paste0(var_name, "_Goterm_UpOrthologs"),   figs_dir, results_dir)
  run_gprofiler(down_ortho$Gene.Group.Identifier,
    paste0(var_name, "_Goterm_DownOrthologs"), figs_dir, results_dir)
  run_gprofiler(c(up_ortho$Gene.Group.Identifier, down_ortho$Gene.Group.Identifier),
    paste0(var_name, "_Goterm_AllOrthologs"),  figs_dir, results_dir)
}

message("\nMLE analysis complete. Output in: ", results_dir)
