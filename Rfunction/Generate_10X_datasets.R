#Generate 10X subsampled data sets

#load 10X data, normalize as in 10X code
wd<-"../10X_full/single-cell-3prime-paper-master/pbmc68k_analysis/"
source(paste(wd,file.path('util.R'),sep=""))
pbmc_68k <- readRDS('../10X_full/data/pbmc68k_data.rds')
all_data <- pbmc_68k$all_data
m<-all_data[[1]]$hg19$mat
l<-.normalize_by_umi(m)  
ExprM.RawCounts.original<-t(l$m)

geneIndices<-l$use_genes
genes<-pbmc_68k$ens_genes[geneIndices]

#load cluster labels
clusters<-read.table(paste(wd,"PureClusters.csv",sep=""),sep=",",header=F) #obtained from 10X code
clusters<-clusters$V1
colnames(ExprM.RawCounts.original)<-clusters
rownames(ExprM.RawCounts.original)=genes

#subsample from data set to create rare clusters
NKcells<-which(clusters=="CD56+ NK")
Bcells<-which(clusters=="CD19+ B")
CD14cells<-which(clusters=="CD14+ Monocyte")

#generate datasets with SAME set of rare cells
#set.seed(50)
samplingRare<-c(sample(NKcells,5))

#subsample to create two common cell types
numbersCD14<-c(800,400,200,100,50,25,13,6)
numbersB<-c(400,200,100,50,25,13,6,3)
for (j in 1:8){
  for (i in 1:20){
    sampling<-c(sample(CD14cells,numbersCD14[j]),sample(Bcells,numbersB[j]))
    sampling<-c(samplingRare,sampling)
    ExprM.RawCounts<-as.matrix(ExprM.RawCounts.original)
    ExprM.RawCounts.sample<-ExprM.RawCounts[,sampling]
    ExprM.RawCounts<-ExprM.RawCounts.sample
    #table(colnames(ExprM.RawCounts.sample))
    
    ExprM.normCounts <- ExprM.RawCounts
    ExprM.RawCounts.filter <- ExprM.RawCounts
    ExprM.normCounts.filter <- ExprM.RawCounts
    
    exprimentID<-paste("10X_rare",j,"_",i,sep="")
    
    #write.table(ExprM.RawCounts, file=paste("results/", exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    #write.table(ExprM.normCounts, file=paste("results/", exprimentID, "_normCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    save(ExprM.RawCounts, file=paste("results/",exprimentID, "_ExprM.RData",sep=""))
    save(ExprM.RawCounts.filter, file=paste("results/", exprimentID, "_ExprM.filter.RData", sep=""))
  }
  
}

