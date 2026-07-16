# sgRNA Analysis Pipeline

## About

Custom pipeline for sgRNA alignment and analysis, developed for cross-disciplinary bioinformatics support at CR2TI, UMR 1064, Nantes Université.

## Prerequisites and genome version
Genome version and annotation :
- sgRNA table used for guide design [sgRNA](sgRNA_Table.tsv)
- [control](Control_sgRNA.tsv) sgRNA table 

## Key analyses
- sgRNA sequence alignment
- On-target / off-target analysis
- Results visualisation

## Prerequisites and genome version
Genome version and annotation :
- sgRNA table used for guide design [sgRNA](sgRNA_Table.tsv)
- [control](Control_sgRNA.tsv) sgRNA table 
- sample annotation for [MLE](Input/SampleAnnot_binary.tsv) or [RRA](Input/SampleAnnot.tsv) models
- Human [orthologs](Input/HumanOrthologs.tsv) from Rabbit Genome 
- MaGeck RRA and MLE output in [Input](Input/) folder

## Tools used for analysis (conda environment)
- R / Bash
- Alignment conda environment => [sgRNAEnv.yml](Scripts/sgRNAEnv.yml)
- R libraries used are shown in all R files

## Scripts used for analysis

- sgRNA Alignment [pipeline](Scripts/sgrna_alignment.sh)
- R scripts for [R_Utility](Scripts/sgrna_utils.R)
- R scripts for [PCA](Scripts/sgrna_pca.R)
- R scripts for [Heatmap](Scripts/sgrna_heatmap.R)
- R scripts for [VolcanoPlot](Scripts/sgrna_volcano.R)
- R scripts for [RankPlot](Scripts/sgrna_rankplot.R)
- R scripts for [MLE](Scripts/sgrna_mle.R)
- R scripts for [GoEnrichment](Scripts/sgrna_go_enrichment.R)

## Author
**Valentin FRANCOIS--CAMPION** — [GitHub](https://github.com/FCValentin)
