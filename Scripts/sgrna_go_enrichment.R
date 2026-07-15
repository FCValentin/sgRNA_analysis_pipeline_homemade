# =============================================================================
# sgrna_go_enrichment.R
# -----------------------------------------------------------------------------
# GO term enrichment for up- and down-regulated genes from MAGeCK output.
# Uses gProfiler2 (g:Profiler) with human ortholog mapping.
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
ORTHOLOG_FILE <- file.path(DATA_DIR, "input/HumanOrthologs.tsv")
PROJECT       <- "Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1"
SCORE         <- 0.01     # FDR threshold
LFC_THRESHOLD <- 1.0      # absolute log2 FC threshold

# =============================================================================
# SETUP
# =============================================================================

figs_dir    <- file.path(DATA_DIR, "Figures/GoTerm", PROJECT)
results_dir <- file.path(DATA_DIR, "results/GoTerm", PROJECT)
create_dirs(figs_dir, results_dir)

# =============================================================================
# LOAD & FILTER
# =============================================================================

message("Loading data...")
genes     <- read_tsv(GENE_SUMMARY)
orthologs <- read_tsv(ORTHOLOG_FILE)

up_genes   <- rownames(genes)[genes$pos.score < SCORE & genes$neg.lfc >  LFC_THRESHOLD]
down_genes <- rownames(genes)[genes$neg.score < SCORE & genes$neg.lfc < -LFC_THRESHOLD]

up_genes   <- up_genes[up_genes %in% orthologs$Symbol]
down_genes <- down_genes[down_genes %in% orthologs$Symbol]

up_ortho   <- orthologs[orthologs$Symbol %in% up_genes, ]
down_ortho <- orthologs[orthologs$Symbol %in% down_genes, ]

write_tsv(up_ortho,   file.path(results_dir, "List_UpOrthologs.tsv"))
write_tsv(down_ortho, file.path(results_dir, "List_DownOrthologs.tsv"))
message("Up: ", nrow(up_ortho), " genes | Down: ", nrow(down_ortho), " genes")

# =============================================================================
# GO ENRICHMENT
# =============================================================================

run_gprofiler(up_ortho$Gene.Group.Identifier,
              "Goterm_UpOrthologs",   figs_dir, results_dir)
run_gprofiler(down_ortho$Gene.Group.Identifier,
              "Goterm_DownOrthologs", figs_dir, results_dir)
run_gprofiler(c(up_ortho$Gene.Group.Identifier,
                down_ortho$Gene.Group.Identifier),
              "Goterm_AllOrthologs",  figs_dir, results_dir)

message("Done. Output in: ", results_dir)
