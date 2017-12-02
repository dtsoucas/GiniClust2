#Generate simulated data using RaceID intestinal data

## load class definition and functions
source("RaceID-master/RaceID_class.R")

## input data- from RaceID
x <- read.csv("RaceID-master/transcript_counts_intestine.xls",sep="\t",header=TRUE)
rownames(x) <- x$GENEID
# prdata: data.frame with transcript counts for all genes (rows) in all cells (columns); with rownames == gene ids; remove ERCC spike-ins 
prdata <- x[grep("ERCC",rownames(x),invert=TRUE),-1]

## RaceID
# initialize SCseq object with transcript counts
sc <- SCseq(prdata)
# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000)
# k-means clustering
sc <- clustexp(sc, clustnr=20,bootnr=50,metric="pearson",do.gap=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000)
# compute t-SNE map
sc <- comptsne(sc,rseed=15555)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75)

#seed used to select cell
seed_means = 10
set.seed(seed_means)
ssc_array = prdata
mean_genes = rowMeans(ssc_array)
#bootstrap sampling
ssc = sample(mean_genes,length(mean_genes),replace = TRUE)
####### Lan code ######
ssc[ssc<0.5] = 0               #filter extremely low expressed genes
names(ssc) = rownames(prdata)  #because length(unique(names(ssc))) = 14870
#######################

num_genes = length(ssc)
size_clusters = c(2000,1000,10,6,4,3)


sim_data = matrix(nrow = nrow(prdata),ncol=sum(size_clusters)) #simulated dataset
rownames(sim_data) = rownames(prdata)
colnames(sim_data) = apply(matrix(1:sum(size_clusters)),1,function(x){paste('C',x,sep='')}) #C1: Cell_1

func_lsize <- sc@background$lsize # calculate the parameter 'size' in function rnbinom
func_lvar  <- sc@background$lvar

#seed used to permute gene
seed_permute = seq(from = 10000, by = 100, length.out = length(size_clusters)-1) 
#seed used to sample cells from each gene distribution
seed_sample = matrix(seq(from = 1000,by = 5,length.out = length(size_clusters)*num_genes),
                     nrow=length(size_clusters)) 

permutedGenes<-c()
for (id_c in 1:length(size_clusters)){
  if(id_c == 1){
    mean_expr = ssc #mean expression values for each gene in cluster 1 (1e-6 avoid generating NA)
    lsize_expr = apply(matrix(mean_expr),1,function(x){func_lsize(x,sc)})
    num_id_c = size_clusters[id_c]#cell number
    sim_data[,1:num_id_c]  = t(apply(matrix(1:num_genes),1,function(i){
      set.seed(seed_sample[id_c,i])
      rnbinom(n=num_id_c,mu=mean_expr[i],size=lsize_expr[i])
    }))
  }else if(id_c>1){
    set.seed(seed_permute[id_c-1])
    id_permuted_genes = sample(which(ssc>=10),size = 100,replace=FALSE)  #permute 100 genes
    permutedGenes<-c(permutedGenes,id_permuted_genes)
    set.seed(seed_permute[id_c-1])
    id_low_genes = sample(which(ssc<10),size = 100,replace=FALSE)
    mean_expr = ssc
    mean_expr[sort(c(id_permuted_genes,id_low_genes))] = mean_expr[c(id_permuted_genes,id_low_genes)] 
    lsize_expr = apply(matrix(mean_expr),1,function(x){func_lsize(x,sc)})
    num_id_c = size_clusters[id_c]#cell number
    sim_data[,(sum(size_clusters[1:(id_c-1)])+1):sum(size_clusters[1:id_c])]= t(apply(matrix(1:num_genes),1,function(i){
      set.seed(seed_sample[id_c,i])
      rnbinom(n=num_id_c,mu=mean_expr[i],size=lsize_expr[i])
    }))
  }
}



sim_data[is.na(sim_data)] <- 0

sim_data_final = cbind(GENEID=rownames(sim_data),sim_data)
dir.create(file.path(workdir, "data"), showWarnings = FALSE) #folder to save data
write.table(sim_data_final,paste("data/Data_",paste(size_clusters,collapse = '_'),".xls",sep = ''),
            row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
