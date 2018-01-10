#Create Fano clustering using k-means

load(paste("results/", exprimentID,"_ExprM.filter.RData", sep=""))  #"ExprM.RawCounts.filter"  "ExprM.normCounts.filter"

#calculate Fano factor, choose top 1000 Fano genes
Fano=apply(ExprM.RawCounts.filter,1,function(x) var(x)/mean(x))
FanoDecreasing<-sort(Fano,decreasing=TRUE)
maxFano=names(FanoDecreasing[1:1000])
m<-ExprM.RawCounts.filter[maxFano,]

#pca on high Fano matrix, use top 50 pcs
pca_Fano<-calcul.pca(t(m),50)

#find k using a gap statistic, if requested
if (gap_statistic==TRUE){
  gapStatistic<-clusGap(pca_Fano$pca,FUN=kmeans,K.max=K.max,iter.max=30)
  k<-with(gapStatistic,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
}

#cluster using k-means
#run k-means 20 times with different seeds, pick seed with lowest error
set.seed(10)
seeds<-sample(1:1000,20,replace=FALSE)
errors<-c()
for (seed in seeds){
  set.seed(seed)
  MaxFanoClusters_pca<-kmeans(pca_Fano$pca,k,iter.max=150,algorithm="MacQueen")
  errors<-c(errors, MaxFanoClusters_pca$tot.withinss)
}
seed<-seeds[which.min(errors)]
set.seed(seed)
MaxFanoClusters_pca<-kmeans(pca_Fano$pca,k,iter.max=150,algorithm="MacQueen")

P_F<-MaxFanoClusters_pca$cluster #Fano clustering result

#save results
save(P_F, file=paste("results/", exprimentID,"_FanoClustering.RData",sep=""))

Fano_pca<-pca_Fano$pca
save(Fano_pca,file=paste("results/", exprimentID,"_Fano_distance.RData",sep=""))

save(maxFano,file=paste("results/", exprimentID,"_maxFanoGenes.RData",sep=""))
