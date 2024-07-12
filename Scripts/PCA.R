### PARAMETERS TO CHANGE ###
directory<-"~/These/Analyses/ProjectGitVanessa" #Enter the path to the directory to use
data_matrix<-"input/sgRNAmatrix.tsv" #path to the matrix file
samples<-"input/SampleAnnot.tsv" #path to the sample file
NbPCA_Axis<-3 # Enter the max number of PCA axis you want to show

# Import some (home) functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")
library(ggplot2)

# Set working directory
setwd(directory)
dir.create("Figures")

# Import proteins background
data <- lire(data_matrix)
sampleAnnot <- lire(samples)

# Check if there are problems on sample names
sampleExprInAnnot<-rn(sampleAnnot)[!rn(sampleAnnot)%in%cn(data)]; if(len(sampleExprInAnnot)>0){warning("Error, these samples don't exist in expression data, please check sample names"); print(sampleExprInAnnot)}
sampleAnnotInExpr<-cn(data)[!cn(data)%in%rn(sampleAnnot)]; if(len(sampleAnnotInExpr)>0){warning("these samples don't exist in sample table, they will be removed from expression data"); print(sampleAnnotInExpr)}
ordered_data <- data[order(row.names(data)),row.names(sampleAnnot)]
print("Sample names match is ok")
exprDat.UPM <- UMI2UPM(ordered_data) #Normalization based on 1M count per sample

##### PCA Quality control#####
print(paste("Total of clusters :",NbPCA_Axis))

Day <- "Day"
FOLR1 <- "FOLR1"
Medium <- "CultureMedium"

acp<-ACP(exprDat.UPM[,]) 
dir.create("Figures/PCA")
pdf(file = paste0("Figures/PCA/PCA_",NbPCA_Axis,"_axis.pdf"),width=10,height=10)
barplot(acp$percentVar,names.arg = round(acp$percentVar*100,2),main = "Contribution of each componant in PCA")

print("Performing PCA....")

for(i in 1:(NbPCA_Axis-1)){
  for(j in (i+1):NbPCA_Axis){
    acp2d(acp,group = as.factor(sampleAnnot[,Day]),plotText = TRUE,pointSize = 4,comp = c(i,j),main="PCA Day",fixedCoord = F)
    acp2d(acp,group = as.factor(sampleAnnot[,Day]),pointSize = 2,comp = c(i,j),main="PCA Day",fixedCoord = F)
    
    acp2d(acp,group = sampleAnnot[,FOLR1],plotText = TRUE,pointSize = 4,comp = c(i,j),main="PCA Naivity",fixedCoord = F)
    acp2d(acp,group = sampleAnnot[,FOLR1],pointSize = 2,comp = c(i,j),main="PCA Naivity",fixedCoord = F)
    
    acp2d(acp,group = sampleAnnot[,Medium],plotText = TRUE,pointSize = 4,comp = c(i,j),main="PCA Medium",fixedCoord = F)
    acp2d(acp,group = sampleAnnot[,Medium],pointSize = 2,comp = c(i,j),main="PCA Medium",fixedCoord = F)
    
    acp2d(acp,pointSize = 2,comp = c(i,j),plotVars = TRUE, plotText = TRUE,fixedCoord = F)
  }
}
dev.off()
print("End of the script, output in Figures folder")



