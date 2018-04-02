#Generate 10X subsampled data sets

#load 10X data, normalize as in 10X code
wd<-"../10X_full/single-cell-3prime-paper-master/pbmc68k_analysis/"
source(paste(wd,file.path('util.R'),sep=""))
pbmc_68k <- readRDS('../10X_full/data/pbmc68k_data.rds')
all_data <- pbmc_68k$all_data
pure_11 <- readRDS('../10X_full/data/all_pure_select_11types.rds')
purified_ref_11 <- load_purified_pbmc_types(pure_11,pbmc_68k$ens_genes)
m<-all_data[[1]]$hg19$mat
l<-.normalize_by_umi(m)  
ExprM.RawCounts.original<-t(l$m)
m_n<-l$m
df<-.get_variable_gene(m_n) 
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off

geneIndices<-l$use_genes
genes<-pbmc_68k$ens_genes[geneIndices]

#load cluster labels
clusters<-read.table(paste(wd,"PureClusters.csv",sep=""),sep=",",header=F) #obtained from 10X code
clusters<-clusters$V1
colnames(ExprM.RawCounts.original)<-clusters
rownames(ExprM.RawCounts.original)=genes

#subsample from data set to create rare clusters

#further refine labels
#NK cells- only consider non-NKT cells- use CD3D, CD3E as NKT cell markers
cd3d<-which(all_data[[1]]$hg19$gene_symbols=="CD3D")
cd3e<-which(all_data[[1]]$hg19$gene_symbols=="CD3E")
cd3g<-which(all_data[[1]]$hg19$gene_symbols=="CD3G")
cd68<-which(all_data[[1]]$hg19$gene_symbols=="CD68")
cd37<-which(all_data[[1]]$hg19$gene_symbols=="CD37")

ind3<-intersect(which(clusters=="CD56+ NK"),which(all_data[[1]]$hg19$mat[,cd3d]==0))
ind4<-intersect(which(clusters=="CD56+ NK"),which(all_data[[1]]$hg19$mat[,cd3e]==0))
ind5<-intersect(which(clusters=="CD56+ NK"),which(all_data[[1]]$hg19$mat[,cd3g]==0))
ind345<-intersect(intersect(ind3,ind4),ind5)


m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]

m_filt<-m_n_1000
use_genes_n<-order(-df$dispersion_norm)
use_genes_n_id<-all_data[[1]]$hg19$gene_symbols[l$use_genes][order(-df$dispersion_norm)]
use_genes_n_ens<-all_data[[1]]$hg19$genes[l$use_genes][order(-df$dispersion_norm)]
z<-.compare_by_cor(m_filt,use_genes_n_ens[1:1000],purified_ref_11)

ind345corrected<-ind345[which(z[ind345,"CD56+ NK"]>quantile(z[ind345,"CD56+ NK"],0.5))]


ind6<-intersect(which(clusters=="CD14+ Monocyte"),which(all_data[[1]]$hg19$mat[,cd68]>1))
ind7<-intersect(which(clusters=="CD14+ Monocyte"),which(all_data[[1]]$hg19$mat[,cd37]>1))
ind67<-intersect(ind6,ind7)
#Monocyte- only consider macrophage cells

NKcells<-ind345corrected
Bcells<-which(clusters=="CD19+ B")
CD14cells<-ind67

#generate datasets with SAME set of rare cells
set.seed(20)
samplingRare<-c(sample(CD14cells,5))
set.seed(100)
seeds<-sample(c(1:10000),140)
#subsample to create two common cell types
numbersNK<-c(1600,800,400,200,100,50,25)
numbersB<-c(800,400,200,100,50,25,13)
for (q in 1:7){
  for (i in 1:20){
    set.seed(seeds[(q-1)*20+i])
    sampling<-c(sample(NKcells,numbersNK[q]),sample(Bcells,numbersB[q]))
    sampling<-c(samplingRare,sampling)
    ExprM.RawCounts<-as.matrix(ExprM.RawCounts.original)
    ExprM.RawCounts.sample<-ExprM.RawCounts[,sampling]
    ExprM.RawCounts<-ExprM.RawCounts.sample
    #table(colnames(ExprM.RawCounts.sample))
    
    ExprM.normCounts <- ExprM.RawCounts
    ExprM.RawCounts.filter <- ExprM.RawCounts
    ExprM.normCounts.filter <- ExprM.RawCounts
    
    exprimentID<-paste("10X_rare",q,"_",i,sep="")
    
    #write.table(ExprM.RawCounts, file=paste("results/", exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    #write.table(ExprM.normCounts, file=paste("results/", exprimentID, "_normCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    save(ExprM.RawCounts, file=paste("results/",exprimentID, "_ExprM.RData",sep=""))
    save(ExprM.RawCounts.filter, file=paste("results/", exprimentID, "_ExprM.filter.RData", sep=""))
  }
  
}

