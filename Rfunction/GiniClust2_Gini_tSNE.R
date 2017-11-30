# Visualization of GiniClust clustering results using tsne 

load(file=paste("results/",exprimentID,"_Gini_distance.RData",sep=""))
load(file=paste("results/", exprimentID,"_GiniClustering.RData",sep=""))

#preliminary pca recommended for larger data sets
if (dim(cell.cell.jaccard.distance)[1]>10000){
  
  pc<-propack.svd(as.matrix(cell.cell.jaccard.distance), neig = 50)
  pcloadings<-t(pc$d*t(pc$u))
  Rtsne_map <- Rtsne(pcloadings, pca = TRUE, perplexity = perplexity_G,check_duplicates = FALSE,max_iter=max_iter_G)  
  
} else{
    
  Rtsne_map <- Rtsne(as.matrix(cell.cell.jaccard.distance), pca = FALSE, perplexity = perplexity_G,check_duplicates = FALSE,max_iter=max_iter_G)       
  
}

ggplot(as.data.frame(Rtsne_map$Y),aes(Rtsne_map$Y[,1],Rtsne_map$Y[,2],col=as.factor(P_G)))+geom_point(size=1,alpha=0.8)+theme_classic()+
  scale_color_manual(values=mycols,name="Gini Clusters")+labs(x="Gini 1", y="Gini 2")
ggsave(paste("figures/Figures_", exprimentID, "_perplexity_", perplexity_G, "_tsne_plot_P_G.pdf", sep=""))
  

# save results
Rtnse_coord2 <-  as.data.frame(Rtsne_map$Y)
rownames(Rtnse_coord2) = rownames(cell.cell.jaccard.distance)
colnames(Rtnse_coord2) = c("dim1", "dim2")
write.table(Rtnse_coord2, file=paste("results/", exprimentID,"_Rtnse_Gini_coord2.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
