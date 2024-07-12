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
library(gprofiler2)

# Set working directory
setwd(directory)
dir.create("results")
dir.create("results/GoTerm")
dir.create(paste0("results/GoTerm/",Project))
dir.create("Figures")
dir.create("Figures/GoTerm")
dir.create(paste0("Figures/GoTerm/",Project))

Genes<-lire(gene_summary)
UpGenes<-row.names(Genes[Genes$pos.score<Score&Genes$neg.lfc>logFCthreshold,])
DownGenes<-row.names(Genes[Genes$neg.score<Score&Genes$neg.lfc<(-logFCthreshold),])

Orthologs<-lire("input/HumanOrthologs.tsv")
UpGenes<-UpGenes[which(UpGenes%in%Orthologs$Symbol)]
DownGenes<-DownGenes[which(DownGenes%in%Orthologs$Symbol)]

UpOrtho<-Orthologs[which(Orthologs$Symbol%in%UpGenes),]
DownOrtho<-Orthologs[which(Orthologs$Symbol%in%DownGenes),]
ecrire(UpOrtho,paste0("results/GoTerm/",Project,"/List_UpOrthologs.tsv"))
ecrire(DownOrtho,paste0("results/GoTerm/",Project,"/List_DownOrthologs.tsv"))

##--------------------------------------------------------##
#### Gprofiler test on R (can also be performed online) ####
##--------------------------------------------------------##

pdf(paste0("Figures/GoTerm/",Project,"/Goterm_UpOrthologs.pdf"),width=10,height=10)
Upgostres <- gost(query = na.omit(UpOrtho$Gene.Group.Identifier),significant = FALSE, organism = "hsapiens")
ecrire(as.matrix(Upgostres$result),paste0("results/GoTerm/",Project,"/Goterm_UpOrthologs.tsv"))
gostplot(Upgostres, capped = FALSE, interactive = FALSE)
dev.off()

pdf(paste0("Figures/GoTerm/",Project,"/Goterm_DownOrthologs.pdf"),width=10,height=10)
Downgostres <- gost(query = na.omit(DownOrtho$Gene.Group.Identifier), significant = FALSE, organism = "hsapiens")
ecrire(as.matrix(Downgostres$result),paste0("results/GoTerm/",Project,"/Goterm_DownOrthologs.tsv"))
gostplot(Downgostres, capped = FALSE, interactive = FALSE)
dev.off()

pdf(paste0("Figures/GoTerm/",Project,"/Goterm_AllOrthologs.pdf"),width=10,height=10)
gostres <- gost(query = na.omit(c(UpOrtho$Gene.Group.Identifier,DownOrtho$Gene.Group.Identifier)), significant = FALSE, organism = "hsapiens")
ecrire(as.matrix(gostres$result),paste0("results/GoTerm/",Project,"/Goterm_AllOrthologs.tsv"))
gostplot(gostres, capped = FALSE, interactive = FALSE)
dev.off()
