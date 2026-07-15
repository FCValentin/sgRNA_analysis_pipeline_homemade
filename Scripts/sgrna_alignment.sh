#!/usr/bin/env bash
# =============================================================================
# sgrna_alignment.sh
# -----------------------------------------------------------------------------
# sgRNA CRISPR screen alignment pipeline using BWA + TrimGalore + MAGeCK.
# Runs as an SGE array job on HPC cluster.
#
# Workflow:
#   1. Build BWA index on sgRNA reference FASTA
#   2. Adapter/quality trimming with TrimGalore + FastQC
#   3. Short-read BWA alignment on sgRNA reference
#   4. Alignment QC (flagstat + idxstats per sample)
#   5. Count sgRNA reads with MAGeCK count
#   6. Differential enrichment tests with MAGeCK test
#   7. MLE modelling with MAGeCK mle
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/sgrna-alignment-pipeline
# Project : sgRNA CRISPR screen — GeCKO library (Vanessa collaboration)
# Date    : 2023 (CR2TI, UMR 1064, Nantes Universite)
# =============================================================================

#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -q max-1m.q
#$ -e ./log/
#$ -o ./log/

set -euo pipefail   # exit on error, undefined variable, or pipe failure

# =============================================================================
# PARAMETERS — edit here
# =============================================================================

REF_FASTA="S014969_sgRNA_rabbit.fasta"
REF_PREFIX="Rabbit_sgRNA"
SGRNA_TABLE="sgRNA_Table.tsv"
SAMPLE_ANNOT="input/SampleAnnot_binary.tsv"
SGRNA_MATRIX="input/sgRNAmatrix.tsv"

TRIM_QUALITY=20
TRIM_LENGTH=20
BWA_THREADS=8
BWA_MIN_SEED=2       # -k: minimum seed length (short for sgRNA ~20 bp)
BWA_GAP_OPEN=0       # -O: gap open penalty

MAGECK_SGRNA_LEN=20
MAGECK_TRIM5=0

# Sample labels for MAGeCK count (comma-separated, matching fastq order)
MAGECK_LABELS="GeCKO_1,GeCKO_2,GeCKO_3,GeCKO_4,GeCKO_5,GeCKO_6,GeCKO_7"

# =============================================================================
# I. DIRECTORY SETUP
# =============================================================================

mkdir -p trimgalore bwa log input

# =============================================================================
# II. BWA INDEX
# =============================================================================

echo "[1/5] Building BWA index on: ${REF_FASTA}"
conda activate VanessaPipeline
bwa index -a is "${REF_FASTA}" -p "${REF_PREFIX}"
conda deactivate

# =============================================================================
# III. TRIMMING (TrimGalore + FastQC)
# =============================================================================

echo "[2/5] Trimming reads..."
conda activate CutRun_ChenPipeline

for sample in $(ls -1 fastq/*_1.fastq.gz | \
               cut -d "/" -f 2 | \
               sed 's/_1.fastq.gz//'); do
    echo "  Trimming: ${sample}"
    trim_galore \
        --fastqc \
        --quality  "${TRIM_QUALITY}" \
        --length   "${TRIM_LENGTH}" \
        -o trimgalore \
        "fastq/${sample}_1.fastq.gz"
done
conda deactivate

# =============================================================================
# IV. BWA ALIGNMENT
# =============================================================================

echo "[3/5] Aligning trimmed reads to sgRNA reference..."
conda activate ChIP_SeqPipeline

for sample in $(ls -1 trimgalore/*_1_trimmed.fq.gz | \
               cut -d "/" -f 2 | \
               sed 's/_1_trimmed.fq.gz//'); do
    echo "  Aligning: ${sample}"
    bwa mem \
        -k "${BWA_MIN_SEED}" \
        -O "${BWA_GAP_OPEN}" \
        -t "${BWA_THREADS}" \
        "${REF_PREFIX}" \
        "trimgalore/${sample}_1_trimmed.fq.gz" \
    | samtools view -hbS \
    | samtools sort \
    > "bwa/${sample}.sort.bam"
done

# =============================================================================
# V. ALIGNMENT QC
# =============================================================================

echo "[4/5] Alignment QC..."

# flagstat per sample
for bam in bwa/*.sort.bam; do
    sample=$(basename "${bam}" .sort.bam)
    samtools flagstat "${bam}" > "log/${sample}.txt"
    echo "  flagstat: ${sample}"
done

# index + idxstats per BAM
for bam in bwa/*.bam; do
    sample=$(basename "${bam}" .bam)
    samtools index  "${bam}"
    samtools idxstats "${bam}" > "log/${sample}.tsv"
    echo "  idxstats: ${sample}"
done

conda deactivate

# =============================================================================
# VI. MAGECK COUNT
# =============================================================================

echo "[5/5] Running MAGeCK count..."
conda activate VanessaPipeline

mageck count \
    -l   "${SGRNA_TABLE}" \
    -n   Test \
    --sgrna-len  "${MAGECK_SGRNA_LEN}" \
    --trim-5     "${MAGECK_TRIM5}" \
    --sample-label "${MAGECK_LABELS}" \
    --fastq \
        trimgalore/cutadapt/cutadapt_GeCKO-31_S1_1.fq.gz \
        trimgalore/cutadapt/cutadapt_GeCKO-32_S2_1.fq.gz \
        trimgalore/cutadapt/cutadapt_GeCKO-33_S3_1.fq.gz \
        trimgalore/cutadapt/cutadapt_GeCKO-34_S4_1.fq.gz \
        trimgalore/cutadapt/cutadapt_GeCKO-35_S5_1.fq.gz \
        trimgalore/cutadapt/cutadapt_GeCKO-36_S6_1.fq.gz \
        trimgalore/cutadapt/cutadapt_GeCKO-37_S7_1.fq.gz

# =============================================================================
# VII. MAGECK TEST (pairwise comparisons)
# =============================================================================

echo "Running MAGeCK test comparisons..."

run_mageck_test() {
    local treatment="$1"
    local control="$2"
    local prefix="$3"
    echo "  Test: ${treatment} vs ${control}"
    mageck test \
        -k "${SGRNA_MATRIX}" \
        -t "${treatment}" \
        -c "${control}" \
        -n "input/${prefix}" \
        --remove-zero both \
        --remove-zero-threshold 0
}

run_mageck_test "D0_None_None"   "D8_Algox_EGFP"   "D0_vs_Day_8_Algox_EGFP"
run_mageck_test "D0_None_None"   "D8_Algox_FOLR1"  "D0_vs_Day_8_Algox_FOLR1"
run_mageck_test "D8_Algox_EGFP"  "D8_Algox_FOLR1"  "Day_8_Algox_EGFP_vs_Day_8_Algox_FOLR1"
run_mageck_test "D8_Algox_EGFP"  "D15_Algox_EGFP"  "Day_8_Algox_EGFP_vs_Day_15_Algox_EGFP"
run_mageck_test "D8_Algox_FOLR1" "D15_Algox_FOLR1" "Day_8_Algox_FOLR1_vs_Day_15_Algox_FOLR1"
run_mageck_test "D8_Algox_EGFP"  "D8_UFO_EGFP"     "Day_8_Algox_EGFP_vs_Day_8_UFO_EGFP"
run_mageck_test "D8_Algox_FOLR1" "D8_UFO_FOLR1"    "Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1"

# =============================================================================
# VIII. MAGECK MLE
# =============================================================================

echo "Running MAGeCK MLE..."
mageck mle \
    --count-table   "${SGRNA_MATRIX}" \
    --design-matrix "${SAMPLE_ANNOT}" \
    --norm-method   total \
    --output-prefix input/test.mle

conda deactivate
echo "Pipeline complete."
