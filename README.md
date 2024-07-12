# Vanessa_sgRNA

### Prerequisites and genome version
Genome version and annotation :
* sgRNA table used for guide design [sgRNA](sgRNA_Table.tsv)
* [control](Control_sgRNA.tsv) sgRNA table 

### Prerequisites and genome version
Genome version and annotation :
* sgRNA table used for guide design [sgRNA](sgRNA_Table.tsv)
* [control](Control_sgRNA.tsv) sgRNA table 
* sample annotation for [MLE](Input/SampleAnnot_binary.tsv) or [RRA](Input/SampleAnnot.tsv) models
* Human [orthologs](Input/HumanOrthologs.tsv) from Rabbit Genome 
* MaGeck RRA and MLE output in [Input](Input/) folder

### Tools used for analysis (conda environment)

* Alignment conda environment => [VanessaEnv.yml](Scripts/VanessaEnv.yml)
* R libraries used are shown in all R files

### Scripts used for analysis

* sgRNA Alignment [pipeline](Scripts/VanessaAlignment.sh)
* R scripts for [PCA](Scripts/PCA.R)
* R scripts for [Heatmap](Scripts/Heatmap.R)
* R scripts for [VolcanoPlot](Scripts/VolcanoPlot.R)
* R scripts for [RankPlot](Scripts/RankPlot.R)
* R scripts for [MLE](Scripts/MLE.R)
* R scripts for [GoEnrichment](Scripts/GoEnrichment.R)

## Authors

** Valentin FRANCOIS--CAMPION / Vanessa QUILLAUD-CHENOUARD ** 
