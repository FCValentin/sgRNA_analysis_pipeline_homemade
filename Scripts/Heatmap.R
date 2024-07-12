### PARAMETERS TO CHANGE ###
directory<-"~/These/Analyses/ProjectGitVanessa" #Enter the path to the directory to use
gene_summary<-"input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.gene_summary.txt" #path to the gene summary file 
sgrna_summary<-"input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.sgrna_summary.txt" #path to the sgRNA summary file
data_matrix<-"input/sgRNAmatrix.tsv" #path to the matrix file
samples<-"input/SampleAnnot.tsv" #path to the sample file
Project<-"Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1" #Name to save the project
reference<-"D0_None_None" #Sample where we need at least one sgRNA expressed as control
Score<-0.01
NbMinSgRNAtoDetect<-2 # 2 or more
SamplesToCompare<-c("D8_Algox_FOLR1","D8_UFO_FOLR1") # Names of sample to compare for manual heatmap

# Import some home functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
library(circlize)
library(ComplexHeatmap)

# Set working directory
setwd(directory)
dir.create("Figures")
dir.create("results")
dir.create("Figures/heatmaps")
dir.create("results/heatmaps")

# Import matrix and samples
data <- lire(data_matrix)
sampleAnnot <- lire(samples)

# Check if there are problems on sample names
sampleExprInAnnot<-rn(sampleAnnot)[!rn(sampleAnnot)%in%cn(data)]; if(len(sampleExprInAnnot)>0){warning("Error, these samples don't exist in expression data, please check sample names"); print(sampleExprInAnnot)}
sampleAnnotInExpr<-cn(data)[!cn(data)%in%rn(sampleAnnot)]; if(len(sampleAnnotInExpr)>0){warning("these samples don't exist in sample table, they will be removed from expression data"); print(sampleAnnotInExpr)}
ordered_data <- data[order(row.names(data)),row.names(sampleAnnot)]
print("Sample names match is ok")
exprDat.UPM <- UMI2UPM(ordered_data) #Normalization based on 1M count per sample

## --------------------------------------------##
#### Script for manual candidate sample list ####
## --------------------------------------------##

Genes<-lire(gene_summary)
Genes_matrix<-exprDat.UPM[data$Gene%in%row.names(Genes)[which((Genes$neg.score<Score|Genes$pos.score<Score)&(Genes$neg.goodsgrna>NbMinSgRNAtoDetect|Genes$pos.goodsgrna>NbMinSgRNAtoDetect))],]
Genes_matrix<-Genes_matrix[which(Genes_matrix[,reference]>0),]

Genes_mat<-(Genes_matrix[,c(SamplesToCompare[1])]+0.001)/(Genes_matrix[,c(SamplesToCompare[2])]+0.001)
Genes_mat<-Genes_mat
quantile.expr<-quantile(unlist(Genes_mat),seq(0,1,.01))
colHA<-colorRamp2(c(0,1,2),c("blue","white","red"))
Ht<-Heatmap(Genes_mat,row_labels = data[row.names(Genes_matrix),"Gene"], row_names_gp = autoGparFontSizeMatrix(nrow(Genes_mat)),
            cluster_rows = FALSE,col = colHA, show_column_names = FALSE,
            cluster_columns = FALSE,name=paste0(SamplesToCompare[1],"/",SamplesToCompare[2]),column_names_gp = autoGparFontSizeMatrix(ncol(Genes_mat)) )

max<-length(Genes_mat)
HeigthScaling<-10
if(max<400){
  HeigthScaling<-10;
}else if(max<1000){
  HeigthScaling<-15;
}else if(max<2000){
  HeigthScaling<-20;
}else if(max<4000){
  HeigthScaling<-40;
}else if(max<6000){
  HeigthScaling<-50;
}else{
  HeigthScaling<-60;
}

dir.create(paste0("Figures/heatmaps/",Project))
dir.create(paste0("results/heatmaps/",Project))

pdf(paste0("Figures/heatmaps/",Project,"/Genes_summary_heatmaps_Score_",Score,".pdf"),width = 10,height=HeigthScaling)
print(Ht)
dev.off()
ecrire(Genes_mat,paste0("results/heatmaps/",Project,"/Genes_summary_heatmaps_Score_",Score,".tsv"))


## -----------------------------------##
#### Script for full sample Heatmap ####
## -----------------------------------##

Genes<-lire(gene_summary)
Genes_matrix<-exprDat.UPM[data$Gene%in%row.names(Genes)[which((Genes$neg.score<Score|Genes$pos.score<Score)&(Genes$neg.goodsgrna>NbMinSgRNAtoDetect|Genes$pos.goodsgrna>NbMinSgRNAtoDetect))],]
Genes_matrix<-Genes_matrix[which(Genes_matrix[,reference]>0),]

Genes_mat<-rowScale(Genes_matrix,center = TRUE,scaled = TRUE)
hclustGeneDE<-unsupervisedClustering(Genes_mat,transpose = F,nboot=30,bootstrap = FALSE,method.dist="pearson")
quantile.expr<-quantile(unlist(Genes_mat),seq(0,1,.01))
colHA<-colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
Ht<-Heatmap(Genes_mat,row_labels = data[row.names(Genes_mat),"Gene"], row_names_gp = autoGparFontSizeMatrix(nrow(Genes_mat)),
            cluster_rows = hclustGeneDE,col = colHA,
            cluster_columns = FALSE,name="Genes_matrix",column_names_gp = autoGparFontSizeMatrix(ncol(Genes_mat)) )

max<-nrow(Genes_mat)
HeigthScaling<-10
if(max<400){
  HeigthScaling<-10;
}else if(max<1000){
  HeigthScaling<-15;
}else if(max<2000){
  HeigthScaling<-20;
}else if(max<4000){
  HeigthScaling<-40;
}else if(max<6000){
  HeigthScaling<-50;
}else{
  HeigthScaling<-60;
}

dir.create(paste0("Figures/heatmaps/",Project))
dir.create(paste0("results/heatmaps/",Project))

pdf(paste0("Figures/heatmaps/",Project,"/Fullsample_Genes_summary_heatmaps_Score_",Score,".pdf"),width = 10,height=HeigthScaling)
print(Ht)
dev.off()
ecrire(Genes_mat,paste0("results/heatmaps/",Project,"/Fullsample_Genes_summary_heatmaps_Score_",Score,".tsv"))

## -----------------------------##
#### Script for sgRNA Heatmap ####
## -----------------------------##

sgrna<-lire(sgrna_summary)
sgrna_matrix<-exprDat.UPM[data$Gene%in%row.names(sgrna)[which(sgrna$p.twosided<Score)],]
sgrna_matrix<-sgrna_matrix[which(sgrna_matrix[,reference]>0),]

sgrna_mat<-rowScale(sgrna_matrix,center = TRUE,scaled = TRUE)
hclustGeneDE<-unsupervisedClustering(sgrna_mat,transpose = F,nboot=30,bootstrap = FALSE,method.dist="pearson")
quantile.expr<-quantile(unlist(sgrna_mat),seq(0,1,.01))
colHA<-colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
Ht<-Heatmap(sgrna_mat,row_labels = data[row.names(sgrna_mat),"Gene"], row_names_gp = autoGparFontSizeMatrix(nrow(sgrna_mat)),
            cluster_rows = hclustGeneDE,col = colHA,
            cluster_columns = FALSE,name="sgrna_matrix",column_names_gp = autoGparFontSizeMatrix(ncol(sgrna_mat)) )

max<-nrow(sgrna_mat)
HeigthScaling<-10
if(max<400){
  HeigthScaling<-10;
}else if(max<1000){
  HeigthScaling<-15;
}else if(max<2000){
  HeigthScaling<-20;
}else if(max<4000){
  HeigthScaling<-40;
}else if(max<6000){
  HeigthScaling<-50;
}else{
  HeigthScaling<-60;
}

pdf(paste0("Figures/heatmaps/",Project,"/Fullsample_sgRNA_summary_heatmaps_Score_",Score,".pdf"),width = 10,height=HeigthScaling)
print(Ht)
dev.off()
ecrire(sgrna_mat,paste0("results/heatmaps/",Project,"/Fullsample_sgRNA_summary_heatmaps_Score_",Score,".tsv"))

