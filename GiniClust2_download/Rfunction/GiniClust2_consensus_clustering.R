#Combine Gini and Fano clustering results using a cluster-aware weighted consensus clustering approach

load(file=paste("results/", exprimentID,"_GiniClustering.RData",sep="")) #P_G
load(file=paste("results/", exprimentID,"_FanoClustering.RData",sep="")) #P_F

P_G_original<-P_G
P_G[which(P_G=="Singleton")]<-NA

#create connectivity matrix from the Gini clustering
M_P_G<-clust2Mat(P_G)
diag(M_P_G)<-1

#to view this heatmap:
#x_hm <- heatmap.2(a, Rowv=FALSE, Colv=FALSE, dendrogram="none", main="8 X 8 Matrix Using Heatmap.2", xlab="Columns", ylab="Rows", col=pal, tracecol="#303030", trace="none", notecol="black", notecex=0.8, keysize = 1.5, margins=c(5, 5))

#create connectivity matrix from the Fano clustering
M_P_F<-clust2Mat(P_F)
diag(M_P_F)<-1

#introduce weighting scheme such that rare cells are more highly weighted in Gini clustering,
#common cells are more highly weighted in Fano clustering
MVWeights<-0.1 # f parameter
clusterSizes<-table(P_G)
clusterFractions<-lapply(clusterSizes,function(x) x/length(P_G))
clusterFractionsList<-rep(NA,length(P_G))
for (i in names(clusterFractions)){
  clusterFractionsList[which(P_G==i)]<-clusterFractions[i]
}
clusterFractionsList<-unlist(unname(clusterFractionsList))
gweights<-sapply(clusterFractionsList,FUN=function(x) 1-plogis(x, scale=0.49*(MinPts/length(P_F)), (MinPts/length(P_F))*4))
gweights[is.na(gweights)]<-0
GWeights<-matrix(rep(NA,length(P_G)**2),ncol=length(P_G),nrow=length(P_G))
for (i in 1:dim(GWeights)[1]){
  for (j in 1:dim(GWeights)[2]){
    if(gweights[i]==gweights[j]){
      GWeights[i,j]<-gweights[i]
    }
    else{
      GWeights[i,j]<-max(gweights[i],gweights[j])
    }
  }
}
#normalize Gweights and MVweights
sumWeights<-GWeights+MVWeights
GWeights<-GWeights/sumWeights
MVWeights<-MVWeights/sumWeights
#find Mtilde
M_P_G_noNA<-M_P_G

for (i in 1:dim(M_P_G)[1]){
  for (j in 1:dim(M_P_G)[1]){
    if (is.na(M_P_G[i,j])){
      M_P_G_noNA[i,j]<-0
    }
  }
}

Mtilde<-M_P_G_noNA*GWeights+M_P_F*MVWeights

save(Mtilde,file=paste("results/", exprimentID,"_Mtilde.RData",sep=""))

#cluster Mtilde to create final clustering
#determine number of final clusters

#final k is number of rare clusters + number of common clusters (if no overlap)
#if rare cluster overlaps with a common cluster (>80%), no extra cluster

#define rare cell group as: rarity < (MinPts/number of cells)
clusterSizesNoSingletons<-clusterSizes[which(names(clusterSizes)!="Singleton")]
rareGroups<-names(clusterSizesNoSingletons)[which(clusterSizesNoSingletons/length(P_F)<(4*MinPts/length(P_F)))]
rareGroupNumber<-length(rareGroups)

#check for overlap with P_F groups
overlap<-0
for (i in 1:rareGroupNumber){
  rareGroupIndices<-which(P_G==rareGroups[i])
  for (commonCluster in unique(P_F)){
    intersection<-length(intersect(rareGroupIndices,which(P_F==commonCluster)))
    overlap<-overlap+as.integer((intersection>0.8*length(rareGroupIndices)) & (intersection>0.75*length(which(P_F==commonCluster))))
  }
}

#k for consensus clustering cannot be larger than the rank of the matrix
#calculate Mtilde matrix rank
Mtilde_rank<-rankMatrix(Mtilde)

#set k for consensus clustering
final_k<-min((k+rareGroupNumber-overlap),Mtilde_rank)

#change to cluster on pca for large matrices- will give the same final result
if (dim(Mtilde)[1]>2**16){
  
  pc<-propack.svd(as.matrix(Mtilde), neig = 100)
  pcloadings<-t(pc$d*t(pc$u))
  Mtilde<-pcloadings
}

set.seed(10)
seeds<-sample(1:1000,20,replace=FALSE)
errors<-c()
for (seed in seeds){
  set.seed(seed)
  consensusClusters<-kmeans(Mtilde,centers=final_k)
  errors<-c(errors, consensusClusters$tot.withinss)
}
seed<-seeds[which.min(errors)]
set.seed(seed)
consensusCluster<-kmeans(Mtilde,centers=final_k)
finalCluster<-consensusCluster$cluster
save(finalCluster,file=paste("results/", exprimentID,"_FinalClustering.RData",sep=""))
