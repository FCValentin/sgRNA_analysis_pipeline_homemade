# =============================================================================
# sgrna_pca.R
# -----------------------------------------------------------------------------
# PCA quality control for sgRNA CRISPR screen count data.
# Plots pairwise 2D PCA panels coloured by Day, FOLR1, and Culture Medium.
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/sgrna-alignment-pipeline
# =============================================================================

source("sgrna_utils.R")
.load_pkg("ggplot2")

# =============================================================================
# PARAMETERS — edit here
# =============================================================================

DATA_DIR     <- "."
DATA_MATRIX  <- file.path(DATA_DIR, "input/sgRNAmatrix.tsv")
SAMPLES_FILE <- file.path(DATA_DIR, "input/SampleAnnot.tsv")
N_AXES       <- 3     # number of PCA axes to display pairwise

# Annotation columns to colour by (must match column names in SampleAnnot.tsv)
COLOR_VARS <- list(
  Day    = "Day",
  FOLR1  = "FOLR1",
  Medium = "CultureMedium"
)

# =============================================================================
# SETUP
# =============================================================================

figs_dir <- file.path(DATA_DIR, "Figures/PCA")
create_dirs(figs_dir)

message("Loading data...")
counts <- read_tsv(DATA_MATRIX)
annot  <- read_tsv(SAMPLES_FILE)
check_sample_names(counts, annot)

ordered_counts <- counts[order(rownames(counts)), rownames(annot)]
upm            <- normalise_upm(ordered_counts)

# =============================================================================
# PCA
# =============================================================================

message("Running PCA (", N_AXES, " axes)...")
pca_res  <- prcomp(t(log1p(as.matrix(upm))), center = TRUE, scale. = FALSE)
pct_var  <- round((pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100, 2)
scores   <- as.data.frame(pca_res$x[, seq_len(N_AXES)])
scores   <- cbind(scores, annot[rownames(scores), , drop = FALSE])

pdf_path <- file.path(figs_dir, paste0("PCA_", N_AXES, "_axes.pdf"))
pdf(pdf_path, width = 10, height = 10)

# Scree plot
barplot(pct_var[seq_len(N_AXES)],
        names.arg = paste0("PC", seq_len(N_AXES), "\n", pct_var[seq_len(N_AXES)], "%"),
        main = "PCA — Variance explained per component",
        ylab = "% Variance", col = "steelblue")

# Pairwise 2D scatter plots
for (i in seq_len(N_AXES - 1)) {
  for (j in seq(i + 1, N_AXES)) {
    pc_x  <- paste0("PC", i)
    pc_y  <- paste0("PC", j)
    x_lab <- paste0(pc_x, " (", pct_var[i], "%)")
    y_lab <- paste0(pc_y, " (", pct_var[j], "%)")

    for (var_label in names(COLOR_VARS)) {
      col_col <- COLOR_VARS[[var_label]]
      if (!col_col %in% colnames(scores)) {
        warning("Column not found in annotation: ", col_col, " — skipping.")
        next
      }
      # With sample labels
      p_text <- ggplot(scores,
                       aes_string(x = pc_x, y = pc_y,
                                  colour = col_col, label = "rownames(scores)")) +
        geom_point(size = 4) +
        geom_text(size = 3, vjust = -0.5, show.legend = FALSE) +
        scale_colour_brewer(palette = "Set1", name = var_label) +
        labs(title = paste("PCA —", var_label, ":", pc_x, "vs", pc_y),
             x = x_lab, y = y_lab) +
        theme_bw(base_size = 11) + coord_fixed()
      print(p_text)

      # Without labels (cleaner overview)
      p_clean <- ggplot(scores,
                        aes_string(x = pc_x, y = pc_y, colour = col_col)) +
        geom_point(size = 2) +
        scale_colour_brewer(palette = "Set1", name = var_label) +
        labs(title = paste("PCA —", var_label, ":", pc_x, "vs", pc_y),
             x = x_lab, y = y_lab) +
        theme_bw(base_size = 11) + coord_fixed()
      print(p_clean)
    }
  }
}
dev.off()
message("PCA saved: ", pdf_path)
