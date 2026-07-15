# =============================================================================
# sgrna_utils.R
# -----------------------------------------------------------------------------
# Shared utility functions for sgRNA CRISPR screen analysis.
# Replaces home functions previously sourced from a private GitLab URL.
# Source this file at the top of every analysis script.
#
# Provides: read_tsv(), write_tsv(), check_sample_names(),
#           normalise_upm(), scale_rows(), build_heatmap_sgrna(),
#           export_heatmap_pdf(), run_gprofiler(), pdf_height_from_nrow()
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/sgrna-alignment-pipeline
# Project : sgRNA CRISPR screen — GeCKO library (Vanessa collaboration)
# Date    : 2023 (CR2TI, UMR 1064, Nantes Universite)
# =============================================================================

.load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' not found. Install with install.packages('", pkg, "')")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

#' Read a tab-separated file with row names
#' @param filepath Character. Path to file.
#' @return data.frame.
read_tsv <- function(filepath) {
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  read.table(filepath, header = TRUE, sep = "\t",
             row.names = 1, stringsAsFactors = FALSE, quote = "")
}

#' Write a data frame to tab-separated file
#' @param x data.frame or matrix.
#' @param filepath Character. Output path.
write_tsv <- function(x, filepath) {
  write.table(x, file = filepath, sep = "\t",
              col.names = NA, row.names = TRUE, quote = FALSE)
  message("Saved: ", filepath)
}

#' Create directory tree (silent if exists)
#' @param ... Character paths to create.
create_dirs <- function(...) {
  for (p in c(...)) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
}

#' Check sample name consistency between expression matrix and annotation
#' @param expr data.frame (genes x samples).
#' @param annot data.frame (samples x variables).
#' @return Invisible NULL. Issues warnings on mismatch.
check_sample_names <- function(expr, annot) {
  miss_expr  <- rownames(annot)[!rownames(annot) %in% colnames(expr)]
  miss_annot <- colnames(expr)[!colnames(expr) %in% rownames(annot)]
  if (length(miss_expr) > 0)
    warning("Samples in annotation but not in matrix: ",
            paste(miss_expr, collapse = ", "))
  if (length(miss_annot) > 0)
    warning("Samples in matrix but not in annotation (will be ignored): ",
            paste(miss_annot, collapse = ", "))
  if (length(miss_expr) == 0 && length(miss_annot) == 0)
    message("Sample names: OK")
  invisible(NULL)
}

#' Normalise count matrix to Units Per Million (UPM)
#' @param counts data.frame or matrix (genes x samples). Raw UMI counts.
#' @return data.frame. UPM-normalised matrix.
normalise_upm <- function(counts) {
  as.data.frame(
    sweep(as.matrix(counts), 2, colSums(counts, na.rm = TRUE) / 1e6, FUN = "/")
  )
}

#' Z-score scale rows of a matrix
#' @param mat Matrix or data.frame (genes x samples).
#' @return Matrix.
scale_rows <- function(mat) {
  t(scale(t(as.matrix(mat)), center = TRUE, scale = TRUE))
}

#' Adaptive PDF height based on number of features
#' @param n Integer. Number of rows.
#' @return Numeric. PDF height in inches.
pdf_height_from_nrow <- function(n) {
  if      (n < 400)  10
  else if (n < 1000) 15
  else if (n < 2000) 20
  else if (n < 4000) 40
  else if (n < 6000) 50
  else               60
}

#' Build a z-score scaled ComplexHeatmap for sgRNA data
#'
#' @param mat          Matrix (sgRNAs x samples). Will be z-score scaled.
#' @param row_labels   Character vector. Row labels (gene names).
#' @param panel_name   Character. Heatmap name (legend label).
#' @param cluster_rows Logical or hclust. Row clustering specification.
#' @param cluster_cols Logical. Column clustering. Default FALSE.
#' @param col_scale    colorRamp2 or NULL. Custom colour scale.
#' @return Heatmap object.
build_heatmap_sgrna <- function(mat, row_labels, panel_name,
                                 cluster_rows = TRUE,
                                 cluster_cols = FALSE,
                                 col_scale    = NULL) {
  .load_pkg("circlize")
  .load_pkg("ComplexHeatmap")

  mat_scaled <- scale_rows(mat)
  if (is.null(col_scale)) {
    q         <- quantile(mat_scaled, probs = c(0.01, 0.99), na.rm = TRUE)
    col_scale <- colorRamp2(c(q[1], 0, q[2]), c("blue", "white", "red"))
  }

  # Pearson-distance row clustering when cluster_rows = TRUE
  if (isTRUE(cluster_rows)) {
    d             <- as.dist(1 - cor(t(mat_scaled),
                                      use   = "pairwise.complete.obs",
                                      method = "pearson"))
    cluster_rows  <- hclust(d, method = "average")
  }

  font_r <- max(2, min(8, 200 / max(nrow(mat_scaled), 1)))
  font_c <- max(4, min(10, 200 / max(ncol(mat_scaled), 1)))

  Heatmap(
    mat_scaled,
    name              = panel_name,
    col               = col_scale,
    row_labels        = row_labels,
    row_names_gp      = gpar(fontsize = font_r),
    column_names_gp   = gpar(fontsize = font_c),
    cluster_rows      = cluster_rows,
    cluster_columns   = cluster_cols,
    show_row_names    = TRUE,
    show_column_names = TRUE
  )
}

#' Export a Heatmap (or list of Heatmaps) to PDF
#' @param ht     Heatmap or list of Heatmaps.
#' @param path   Character. PDF output path.
#' @param n_rows Integer. Number of rows (used for adaptive height).
export_heatmap_pdf <- function(ht, path, n_rows) {
  h <- pdf_height_from_nrow(n_rows)
  pdf(path, width = 10, height = h)
  print(ht)
  dev.off()
  message("Heatmap saved: ", path)
}

#' Run gProfiler2 enrichment and save results
#'
#' @param gene_ids   Character vector. Gene identifiers (Ensembl or HGNC).
#' @param label      Character. Label used in output filenames.
#' @param out_figs   Character. Figures directory.
#' @param out_res    Character. Results directory.
#' @param organism   Character. gProfiler organism code. Default "hsapiens".
#' @return Invisible gost result object, or NULL if input is empty.
run_gprofiler <- function(gene_ids, label, out_figs, out_res,
                           organism = "hsapiens") {
  .load_pkg("gprofiler2")
  gene_ids <- na.omit(gene_ids)
  if (length(gene_ids) == 0) {
    message("  No genes for enrichment: ", label, " — skipped.")
    return(invisible(NULL))
  }
  res <- gost(query = gene_ids, significant = FALSE, organism = organism)
  if (!is.null(res)) {
    write_tsv(as.data.frame(res$result),
              file.path(out_res, paste0(label, ".tsv")))
    pdf(file.path(out_figs, paste0(label, ".pdf")), width = 10, height = 10)
    print(gostplot(res, capped = FALSE, interactive = FALSE))
    dev.off()
    message("  GO enrichment saved: ", label)
  }
  invisible(res)
}
