### PARAMETERS TO CHANGE ###
directory<-"~/These/Analyses/ProjectGitVanessa" #Enter the path to the directory to use
data_matrix<-"input/sgRNAmatrix.tsv" #path to the matrix file
gene_summary<-"input/test.mle.gene_summary.txt" #path to the gene summary file 
sgrna_summary<-"input/test.mle.sgrna_summary.txt" #path to the sgRNA summary file
Project<-"MLE" #Name to save the project
samples<-"input/SampleAnnot_binary.tsv" #path to the sample file

reference<-"D0_None_None" #Sample where we need at least one sgRNA expressed as control
Score<-2 #Beta score to test
NbMinSgRNAtoDetect<-2 # 2 or more
SamplesToCompare<-c("D8_Algox_FOLR1","D8_UFO_FOLR1") # Names of sample to compare for manual heatmap

# Import some home functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
library(circlize)
library(ComplexHeatmap)
library(gprofiler2)

# Set working directory
setwd(directory)
dir.create("Figures")
dir.create("results")
dir.create("Figures/MLE")
dir.create("results/MLE")

# Import matrix and samples
data <- lire(data_matrix)
sampleAnnot <- lire(samples)
Genes<-lire(gene_summary)
sgrna<-read.table(sgrna_summary,row.names = 2,header = T)
Orthologs<-lire("input/HumanOrthologs.tsv")

# Check if there are problems on sample names
sampleExprInAnnot<-rn(sampleAnnot)[!rn(sampleAnnot)%in%cn(data)]; if(len(sampleExprInAnnot)>0){warning("Error, these samples don't exist in expression data, please check sample names"); print(sampleExprInAnnot)}
sampleAnnotInExpr<-cn(data)[!cn(data)%in%rn(sampleAnnot)]; if(len(sampleAnnotInExpr)>0){warning("these samples don't exist in sample table, they will be removed from expression data"); print(sampleAnnotInExpr)}
ordered_data <- data[order(row.names(data)),row.names(sampleAnnot)]
print("Sample names match is ok")
exprDat.UPM <- UMI2UPM(ordered_data) #Normalization based on 1M count per sample


for(i in 2:ncol(sampleAnnot)){
  VariableTested<-paste0(colnames(sampleAnnot)[i],".beta") # extract beta score
  Genes_matrix<-exprDat.UPM[data$Gene%in%row.names(Genes)[which((Genes[,VariableTested]<(-Score)|Genes[,VariableTested]>Score))],] # extract genes significant for beta score parameter
  Genes_matrix<-Genes_matrix[which(Genes_matrix[,reference]>0),] # filter samples with no reads in the reference sample
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
  
  pdf(paste0("Figures/MLE/",colnames(sampleAnnot)[i],"_Fullsample_Genes_summary_heatmaps_BetaScore_",Score,".pdf"),width = 10,height=HeigthScaling)
  hist(Genes[,VariableTested],xlab = "BetaScore",main = paste0("Variable tested : ",colnames(sampleAnnot)[i])) # histogram distribution
  print(Ht)
  dev.off()
  ecrire(Genes_mat,paste0("results/MLE/",colnames(sampleAnnot)[i],"_Fullsample_Genes_summary_heatmaps_BetaScore_",Score,".tsv"))
 
  ##-------------##
  #### Go term ####
  ##-------------##
  UpGenes<-row.names(Genes)[which(Genes[,VariableTested]>Score)]
  DownGenes<-row.names(Genes)[which(Genes[,VariableTested]<(-Score))]
  
  UpGenes<-UpGenes[which(UpGenes%in%Orthologs$Symbol)]
  DownGenes<-DownGenes[which(DownGenes%in%Orthologs$Symbol)]
  
  UpOrtho<-Orthologs[which(Orthologs$Symbol%in%UpGenes),]
  DownOrtho<-Orthologs[which(Orthologs$Symbol%in%DownGenes),]
  ecrire(UpOrtho,paste0("results/MLE/",colnames(sampleAnnot)[i],"_List_UpOrthologs.tsv"))
  ecrire(DownOrtho,paste0("results/MLE/",colnames(sampleAnnot)[i],"_List_DownOrthologs.tsv"))
  
  ##--------------------------------------------------------##
  #### Gprofiler test on R (can also be performed online) ####
  ##--------------------------------------------------------##
  
  pdf(paste0("Figures/MLE/",colnames(sampleAnnot)[i],"_Goterm_UpOrthologs.pdf"),width=10,height=10)
  Upgostres <- gost(query = na.omit(UpOrtho$Gene.Group.Identifier),significant = FALSE, organism = "hsapiens")
  ecrire(as.matrix(Upgostres$result),paste0("results/MLE/",colnames(sampleAnnot)[i],"_Goterm_UpOrthologs.tsv"))
  gostplot(Upgostres, capped = FALSE, interactive = FALSE)
  dev.off()
  
  pdf(paste0("Figures/MLE/",colnames(sampleAnnot)[i],"_Goterm_DownOrthologs.pdf"),width=10,height=10)
  Downgostres <- gost(query = na.omit(DownOrtho$Gene.Group.Identifier), significant = FALSE, organism = "hsapiens")
  ecrire(as.matrix(Downgostres$result),paste0("results/MLE/",colnames(sampleAnnot)[i],"_Goterm_DownOrthologs.tsv"))
  gostplot(Downgostres, capped = FALSE, interactive = FALSE)
  dev.off()
  
  pdf(paste0("Figures/MLE/",colnames(sampleAnnot)[i],"_Goterm_AllOrthologs.pdf"),width=10,height=10)
  gostres <- gost(query = na.omit(c(UpOrtho$Gene.Group.Identifier,DownOrtho$Gene.Group.Identifier)), significant = FALSE, organism = "hsapiens")
  ecrire(as.matrix(gostres$result),paste0("results/MLE/",colnames(sampleAnnot)[i],"_Goterm_AllOrthologs.tsv"))
  gostplot(gostres, capped = FALSE, interactive = FALSE)
  dev.off()
}
