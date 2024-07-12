### PARAMETERS TO CHANGE ###
directory<-"~/These/Analyses/ProjectGitVanessa" #Enter the path to the directory to use
gene_summary<-"input/Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1.gene_summary.txt" #path to the gene summary file 
Project<-"Day_8_Algox_FOLR1_vs_Day_8_UFO_FOLR1" #Name to save the project
Score<-0.01
TopRank<-10 #Top rank to take

# Import some home functions
source("https://gitlab.univ-nantes.fr/E114424Z/veneR/raw/master/loadFun.R?inline=false")

# Set working directory
setwd(directory)
dir.create("Figures")
dir.create("Figures/RankPlot")
dir.create(paste0("Figures/RankPlot/",Project))

## ----------------------------------------------------##
#### Script for Rank plot (Neg and then Pos selected ####
## ----------------------------------------------------##

NegGenes<-lire(gene_summary)
pdf(paste0("Figures/RankPlot/",Project,"/Genes_summary_Score_",Score,".pdf"),width=10,height=10)
plot(NegGenes$neg.rank,-log10(NegGenes$neg.score),type="n",ylab="-log10(neg.score)",xlab="rank",col="black", main=paste0("Neg Score ranking : ",Project))
if(nrow(NegGenes[NegGenes$neg.score<Score&NegGenes$neg.rank<TopRank,])>0) text(NegGenes[NegGenes$neg.score<Score&NegGenes$neg.rank<TopRank,]$neg.rank,-log10(NegGenes[NegGenes$neg.score<Score&NegGenes$neg.rank<TopRank,]$neg.score),labels = row.names(NegGenes)[NegGenes$neg.score<Score&NegGenes$neg.rank<TopRank],cex=0.5,col="green")
points(NegGenes$neg.rank,-log10(NegGenes$neg.score),cex=.2,col="black",pch=16)

PosGenes<-NegGenes[order(NegGenes$pos.rank),]
plot(PosGenes$pos.rank,-log10(PosGenes$pos.score),type="n",ylab="-log10(pos.score)",xlab="rank",col="black", main=paste0("Pos Score ranking : ",Project))
if(nrow(PosGenes[PosGenes$pos.score<Score&PosGenes$pos.rank<TopRank,])>0) text(PosGenes[PosGenes$pos.score<Score&PosGenes$pos.rank<TopRank,]$pos.rank,-log10(PosGenes[PosGenes$pos.score<Score&PosGenes$pos.rank<TopRank,]$pos.score),labels = row.names(PosGenes)[PosGenes$pos.score<Score&PosGenes$pos.rank<TopRank],cex=0.5,col="red")
points(PosGenes$pos.rank,-log10(PosGenes$pos.score),cex=.2,col="black",pch=16)
dev.off()
