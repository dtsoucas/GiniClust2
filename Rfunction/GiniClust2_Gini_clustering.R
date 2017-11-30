#Clustering using high Gini genes

#########################################################
### G) binarization, jaccard distance, dbscan for clustering
#########################################################

# .1 load data
load(paste("results/", exprimentID,"_StatisticsTable_afterLOESS.RData",sep=""))  #"Genelist.top_pvalue" "ExprM.Stat2" 
load(paste("results/", exprimentID,"_ExprM.filter.RData", sep=""))  #"ExprM.RawCounts.filter"  "ExprM.normCounts.filter"           
load(paste("results/", exprimentID,"_ExprM.RData", sep=""))  #"ExprM.RawCounts"  "ExprM.normCounts"           

dim(ExprM.RawCounts.filter)

if(GeneList == 1)
{
  GeneList.final       = Genelist.top_pvalue
}else{
  GeneList.final       = Genelist.HighNormGini
}

# .2 binarization
m = ExprM.RawCounts.filter[GeneList.final,]
#m<-as.matrix(m)
m2 = m
# if t is the expression vector, the question asked here is which count value x is the smallest one, when sum(t[t>x]) > Gamma % * sum(t). 
bc.list.low = c()
bc.list.med = c()
bc.list.high = c()
for(rn in 1:dim(m)[1])
{
  t <- as.numeric(m[rn,])
  t.table <- data.frame(table(t))
  c=as.numeric(as.character(t.table$t))
  f=as.numeric(as.character(t.table$Freq))
  tcf  <- data.frame(cbind(c,f))
  csum <- apply(tcf,1,function(x){sum(t[t>=x[1]])/sum(t[t>0])})
  tcfs <- data.frame(cbind(c,f,csum))
  tcfs <- tcfs[rev(order(tcfs$c)),]
  n = max(3, which(tcfs$csum>Gamma)[1])
  bc.list.high[rn]  = tcfs$c[n]
  bc.list.low[rn] = tcfs$c[n+1]   #range will be > bc.list.low[rn] to bc.list.high[rn] 
  bc.list.med[rn] = mean(c(bc.list.high[rn], bc.list.low[rn]))
}
top_n_gene = max(length(bc.list.low)*0.1, 10)
RawCounts_cutoff = floor(mean(bc.list.med[1:top_n_gene]))
all.gene.as.col.ExprM.RawCounts.binary = t(ExprM.RawCounts.filter)
all.gene.as.col.ExprM.RawCounts.binary[all.gene.as.col.ExprM.RawCounts.binary <  RawCounts_cutoff] = 0
all.gene.as.col.ExprM.RawCounts.binary[all.gene.as.col.ExprM.RawCounts.binary >= RawCounts_cutoff] = 1
final.gene.as.col.ExprM.RawCounts.binary = all.gene.as.col.ExprM.RawCounts.binary[, intersect(colnames(all.gene.as.col.ExprM.RawCounts.binary), GeneList.final)]
dim(final.gene.as.col.ExprM.RawCounts.binary) 

# .3 jaccard distance matrix
###------------------------------- Complementary codes -------------------------------###
#locate cells whose gene expression values are all zeros
index.cell.zero = which(apply(final.gene.as.col.ExprM.RawCounts.binary,1,function(x) length(which(x>0)))==0)
if (dim(ExprM.RawCounts.filter)[2]<2**16){
  cell.cell.jaccard.distance = 1 - jaccard(final.gene.as.col.ExprM.RawCounts.binary)
  cell.cell.jaccard.distance = as.data.frame(as.matrix(cell.cell.jaccard.distance))
  #covert distance between two 'zero cells' to zero
  cell.cell.jaccard.distance[index.cell.zero,index.cell.zero] = 0
  cell.cell.jaccard.distance = as.data.frame(as.matrix(cell.cell.jaccard.distance))
  rownames(cell.cell.jaccard.distance)  = make.unique(rownames(final.gene.as.col.ExprM.RawCounts.binary))
  colnames(cell.cell.jaccard.distance)  = rownames(final.gene.as.col.ExprM.RawCounts.binary)
}else{
  cell.cell.jaccard.distance = jaccard_dist_large_matrix(final.gene.as.col.ExprM.RawCounts.binary)
}

#rownames(cell.cell.jaccard.distance)  = rownames(final.gene.as.col.ExprM.RawCounts.binary)
dim(cell.cell.jaccard.distance)
###-----------------------------------------------------------------------------------###

#determine MinPts
if (automatic_minpts==TRUE){
  if (dim(cell.cell.jaccard.distance)[1]>3000){
    #set to 0.1% of total number of points
    MinPts<-floor(.001*dim(cell.cell.jaccard.distance)[1])
  } else {MinPts<-3}
}
#determine eps
if (automatic_eps==TRUE){
  if (dim(cell.cell.jaccard.distance)[1]<10000){ #after 10000 this becomes intensive to calculate
    knndist<-dbscan::kNNdist(cell.cell.jaccard.distance, k =  MinPts)
    knndistNoZero<-knndist
    if (length(which(knndist==0))>0){
      knndistNoZero<-knndist[-which(knndist==0)]
    }
    eps=(sort(knndistNoZero)[floor(.00125*(dim(cell.cell.jaccard.distance)[1])*(MinPts))]) #2057
  } else {
    #compute eps for multiple smaller subsamples, use average eps
    epss<-c()
    for (i in 1:10){
      size<-dim(cell.cell.jaccard.distance)[1]*3/MinPts
      sampleIndices<-sample(seq(1:dim(cell.cell.jaccard.distance)[1]),size,replace=FALSE)
      sampleDistance<-cell.cell.jaccard.distance[sampleIndices,sampleIndices]
      knndist<-dbscan::kNNdist(sampleDistance, k =  3)
      knndistNoZero<-knndist
      if (length(which(knndist==0))>0){
        knndistNoZero<-knndist[-which(knndist==0)]
      }
      eps=mean(sort(knndistNoZero)[1:floor(.00125*(size*(3)))]) #2057
      epss<-c(epss,eps)
    }
    eps = mean(epss)
  }
}

# .3 dbscan
title             = paste("eps", eps, "MinPts", MinPts, sep=".")
data.mclust       = fpc::dbscan(cell.cell.jaccard.distance, eps = eps, MinPts = MinPts,  method = "dist", showplot = FALSE)  

#rename cluster names to include singletons
o_membership      = factor(paste("db_",data.mclust$cluster ,sep=""))
if(levels(o_membership)[1]=="db_0"){levels(o_membership)[1] = "Singleton"}
c_membership = o_membership

P_G=c_membership

save(P_G, file=paste("results/", exprimentID,"_GiniClustering.RData",sep=""))
save(cell.cell.jaccard.distance,file=paste("results/",exprimentID,"_Gini_distance.RData",sep=""))
save(GeneList.final,file=paste("results/",exprimentID,"_GeneList_final.RData",sep=""))
