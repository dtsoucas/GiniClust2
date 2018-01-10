# Visualization of GiniClust clustering results using tsne 

load(file=paste("results/",exprimentID,"_Gini_distance.RData",sep=""))
load(file=paste("results/", exprimentID,"_GiniClustering.RData",sep=""))

#preliminary pca recommended for larger data sets
if (dim(cell.cell.jaccard.distance)[1]>10000){
  
  set.seed(10)
  pc<-propack.svd(as.matrix(cell.cell.jaccard.distance), neig = 50)
  pcloadings<-t(pc$d*t(pc$u))
  Rtsne_map <- Rtsne(pcloadings, pca = TRUE, perplexity = perplexity_G,check_duplicates = FALSE,max_iter=max_iter_G,dims=1)  
  
} else{
    
  set.seed(10)
  Rtsne_map <- Rtsne(as.matrix(cell.cell.jaccard.distance), pca = FALSE, perplexity = perplexity_G,check_duplicates = FALSE,max_iter=max_iter_G,dims=1)       
  
}

# save results
Rtnse_coord2 <-  as.data.frame(Rtsne_map$Y)
rownames(Rtnse_coord2) = rownames(cell.cell.jaccard.distance)
colnames(Rtnse_coord2) = c("dim1")
write.table(Rtnse_coord2, file=paste("results/", exprimentID,"_Rtnse_Gini_coord2.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
