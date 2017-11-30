#Preprocess and filter 10X data using 10X code: https://github.com/10XGenomics/single-cell-3prime-paper

wd<-"single-cell-3prime-paper-master/pbmc68k_analysis/"
source(paste(wd,file.path('util.R'),sep=""))
pbmc_68k <- readRDS(paste("data/",file.path('pbmc68k_data.rds'),sep=""))
all_data <- pbmc_68k$all_data
m<-all_data[[1]]$hg19$mat
genes<-all_data[[1]]$hg19$genes

ExprM.RawCounts<-as.matrix(t(m))
rownames(ExprM.RawCounts)<-genes

l<-.normalize_by_umi(m) 
genes_used<-genes[l$use_genes]


ExprM.RawCounts.filter<-as.matrix(t(l$m))
rownames(ExprM.RawCounts.filter)<-genes_used

save(ExprM.RawCounts, file=paste("results/",exprimentID, "_ExprM.RData",sep=""))
save(ExprM.RawCounts.filter, file=paste("results/", exprimentID, "_ExprM.filter.RData", sep=""))
