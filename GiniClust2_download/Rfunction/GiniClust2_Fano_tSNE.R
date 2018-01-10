# Visualization of Fano clustering using tsne 

load(file=paste("results/",exprimentID,"_Fano_distance.RData",sep=""))
load(file=paste("results/", exprimentID,"_FanoClustering.RData",sep=""))

n<-length(unique(P_F))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


set.seed(10)
Rtsne_map <- Rtsne(as.matrix(Fano_pca), pca = TRUE, perplexity = perplexity_F,check_duplicates = FALSE,max_iter=max_iter_F)       

if (dim(Fano_pca)[1]>1000){
  #plot with smaller points
  
  ggplot(as.data.frame(Rtsne_map$Y),aes(Rtsne_map$Y[,1],Rtsne_map$Y[,2],col=as.factor(P_F)))+geom_point(size=0.5,alpha=0.6)+theme_classic()+
    scale_color_manual(values=mycols,name="Fano Clusters")+labs(x="Fano 1", y="Fano 2")
  ggsave(paste("figures/Figures_", exprimentID, "_perplexity_", perplexity_F, "_tsne_plot_P_F.pdf", sep=""))
} else{
  #plot with larger points
  
  ggplot(as.data.frame(Rtsne_map$Y),aes(Rtsne_map$Y[,1],Rtsne_map$Y[,2],col=as.factor(P_F)))+geom_point(size=1,alpha=0.8)+theme_classic()+
    scale_color_manual(values=mycols,name="Fano Clusters")+labs(x="Fano 1", y="Fano 2")
  
  ggsave(paste("figures/Figures_", exprimentID, "_perplexity_", perplexity_F, "_tsne_plot_P_F.pdf", sep=""))
  
}


# Save results
Rtnse_coord2 <-  as.data.frame(Rtsne_map$Y)
colnames(Rtnse_coord2) = c("dim1", "dim2")
write.table(Rtnse_coord2, file=paste("results/", exprimentID,"_Rtnse_Fano_coord2.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
