### PARAMETERS TO CHANGE ###
directory<-"~/These/Analyses/ProjectGitVanessa" #Enter the path to the directory to use
gene_summary<-"input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.gene_summary.txt" #path to the gene summary file 
sgrna_summary<-"input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.sgrna_summary.txt" #path to the sgRNA summary file
Project<-"Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1" #Name to save the project

data_matrix<-"input/sgRNAmatrix.tsv" #path to the matrix file
Score<-0.01
logFCthreshold<-1 #Absolute Log2(Fold-Change) threshold (if logFCthreshold=1, gene is differentially expressed if expressed 2 time more or less between folds)

# Import some home functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")

# Set working directory
setwd(directory)
dir.create("Figures")
dir.create("Figures/Volcano")
dir.create(paste0("Figures/Volcano/",Project))


## -----------------------------##
#### Script for genes Volcano ####
## -----------------------------##

Genes<-lire(gene_summary)

pdf(paste0("Figures/Volcano/",Project,"/VolcanoPlot_Genes_summary_Score_",Score,".pdf"),width=10,height=10)
plot(Genes$neg.lfc,-log10(Genes$neg.score),type="n",ylab="-log10(neg.score)",xlab="neg.lfc",col="black",
     main=paste0("Volcano plot of comparison : ",Project))
abline(h=-log10(Score))
abline(v = -logFCthreshold)
abline(v = logFCthreshold)
if(nrow(Genes[Genes$pos.score<Score&Genes$neg.lfc>logFCthreshold,])>0) text(Genes[Genes$pos.score<Score&Genes$neg.lfc>logFCthreshold,]$neg.lfc,-log10(Genes[Genes$pos.score<Score&Genes$neg.lfc>logFCthreshold,]$pos.score),labels = row.names(Genes)[Genes$pos.score<Score&Genes$neg.lfc>logFCthreshold],cex=.2,col="red")
if(nrow(Genes[Genes$neg.score<Score&Genes$neg.lfc<(-logFCthreshold),])>0) text(Genes[Genes$neg.score<Score&Genes$neg.lfc<(-logFCthreshold),]$neg.lfc,-log10(Genes[Genes$neg.score<Score&Genes$neg.lfc<(-logFCthreshold),]$neg.score),labels = row.names(Genes)[Genes$neg.score<Score&Genes$neg.lfc<(-logFCthreshold)],cex=0.5,col="green")
points(c(Genes$neg.lfc,Genes$pos.lfc),c(-log10(Genes$neg.score),-log10(Genes$pos.score)),cex=.2,col="black",pch=16)
dev.off()

## -----------------------------##
#### Script for sgRNA Volcano ####
## -----------------------------##

sgrna<-lire(sgrna_summary)
data <- lire(data_matrix)


pdf(paste0("Figures/Volcano/",Project,"/VolcanoPlot_sgRNA_summary_Score_",Score,".pdf"),width=10,height=10)
plot(sgrna$LFC,-log10(sgrna$FDR),type="n",ylab="FDR",xlab="LFC",col="black",
     main=paste0("Volcano plot of comparison : ",Project))
abline(h= -log10(Score))
abline(v = -logFCthreshold)
abline(v = logFCthreshold)
if(nrow(sgrna[sgrna$FDR<Score&sgrna$LFC>logFCthreshold,])>0) text(sgrna[sgrna$FDR<Score&sgrna$LFC>logFCthreshold,]$LFC,-log10(sgrna[sgrna$FDR<Score&sgrna$LFC>logFCthreshold,]$FDR),labels = data[row.names(sgrna)[sgrna$FDR<Score&sgrna$LFC>logFCthreshold],"Gene"],cex=0.3,col="red")
if(nrow(sgrna[sgrna$FDR<Score&sgrna$LFC<(-logFCthreshold),])>0) text(sgrna[sgrna$FDR<Score&sgrna$LFC<(-logFCthreshold),]$LFC,-log10(sgrna[sgrna$FDR<Score&sgrna$LFC<(-logFCthreshold),]$FDR),labels = data[row.names(sgrna)[sgrna$FDR<Score&sgrna$LFC<(-logFCthreshold)],"Gene"],cex=0.3,col="green")
points(sgrna$LFC,-log10(sgrna$FDR),cex=.2,col="black",pch=16)
dev.off()

