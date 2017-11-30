#composite tSNE visualization

#load Gini and Fano dimensions
Rtsne_map_G<-read.table(paste("results/", exprimentID,"_Rtnse_Gini_coord2.csv", sep=""), sep=",")
Rtsne_map_F<-read.table(paste("results/", exprimentID,"_Rtnse_Fano_coord2.csv", sep=""), sep=",")

#load final clustering ids
load(file=paste("results/", exprimentID,"_FinalClustering.RData",sep=""))

#load filtered data set
load(paste("results/", exprimentID,"_ExprM.filter.RData", sep=""))  #"ExprM.RawCounts.filter"  "ExprM.normCounts.filter"


#plot visualizations

if (length(finalCluster)>1000){
  #smaller point size
  p1<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,1],Rtsne_map_G[,1],col=as.factor(finalCluster)))+geom_point(size=0,alpha=0.6)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 1", y="Gini 1")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p2<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,1],Rtsne_map_G[,2],col=as.factor(finalCluster)))+geom_point(size=0,alpha=0.6)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 1", y="Gini 2")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p3<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,2],Rtsne_map_G[,1],col=as.factor(finalCluster)))+geom_point(size=0,alpha=0.6)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 2", y="Gini 1")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p4<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,2],Rtsne_map_G[,2],col=as.factor(finalCluster)))+geom_point(size=0,alpha=0.6)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 2", y="Gini 2")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p5<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,1],Rtsne_map_F[,2],col=as.factor(finalCluster)))+geom_point(size=0,alpha=0.6)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 1", y="Fano 2")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p6<-ggplot(Rtsne_map_F,aes(Rtsne_map_G[,1],Rtsne_map_G[,2],col=as.factor(finalCluster)))+geom_point(size=0,alpha=0.6)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Gini 1", y="Gini 2")+ guides(colour = guide_legend(override.aes = list(size=2)))
  g<-grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
  ggsave(paste("figures/Figures_", exprimentID, "_tsne_plot_combined.pdf", sep=""),g)
  
} else{
  #larger point size
  p1<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,1],Rtsne_map_G[,1],col=as.factor(finalCluster)))+geom_point(size=1,alpha=0.8)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 1", y="Gini 1")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p2<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,1],Rtsne_map_G[,2],col=as.factor(finalCluster)))+geom_point(size=1,alpha=0.8)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 1", y="Gini 2")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p3<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,2],Rtsne_map_G[,1],col=as.factor(finalCluster)))+geom_point(size=1,alpha=0.8)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 2", y="Gini 1")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p4<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,2],Rtsne_map_G[,2],col=as.factor(finalCluster)))+geom_point(size=1,alpha=0.8)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 2", y="Gini 2")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p5<-ggplot(Rtsne_map_F,aes(Rtsne_map_F[,1],Rtsne_map_F[,2],col=as.factor(finalCluster)))+geom_point(size=1,alpha=0.8)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Fano 1", y="Fano 2")+ guides(colour = guide_legend(override.aes = list(size=2)))
  p6<-ggplot(Rtsne_map_F,aes(Rtsne_map_G[,1],Rtsne_map_G[,2],col=as.factor(finalCluster)))+geom_point(size=1,alpha=0.8)+theme_classic()+
    scale_color_manual(values=mycols,name="Final Clusters")+labs(x="Gini 1", y="Gini 2")+ guides(colour = guide_legend(override.aes = list(size=2)))
  g<-grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2)
  ggsave(paste("figures/Figures_", exprimentID, "_tsne_plot_combined.pdf", sep=""),g)
  
}


#venn diagram figures to compare overlap with DE genes and high Fano and high Gini genes
#########################################################

load(file=paste("results/",exprimentID,"_GeneList_final.RData",sep=""))
load(file=paste("results/",exprimentID,"_maxFanoGenes.RData",sep=""))

clusterSizes<-table(finalCluster)
rareGroups<-names(clusterSizes)[which(clusterSizes/length(finalCluster)<(4*MinPts/length(finalCluster)))]
commonGroups<-names(clusterSizes)[which(clusterSizes/length(finalCluster)>=(4*MinPts/length(finalCluster)))]

#for rare groups, compare to Gini
for(rare.cluster in rareGroups ){
  # .1 vennplot plot 
  load(paste("results/",rare.cluster,"_MIST_Sig.RData", sep=""))
  
  # .2 overlap fisher.test
  rownames(cluster_lrTest.table) = cluster_lrTest.table$Gene
  cluster_lrTest.table$log2fold_change = as.numeric(as.character(cluster_lrTest.table$log2fold_change)) 
  cluster.diffgene <- rownames(cluster_lrTest.table[cluster_lrTest.table$p_value < lr.p_value_cutoff  & abs(cluster_lrTest.table$log2fold_change) > diff.cutoff,])
  length(cluster.diffgene)
  length(intersect(cluster.diffgene, GeneList.final)) 
  intersect(cluster.diffgene, GeneList.final)
  
  area1 = length(cluster.diffgene)
  area2 = length(GeneList.final)
  n12   = length(intersect(cluster.diffgene, GeneList.final))
  cluster.overlap <- matrix(c(n12, area1-n12 , area2-n12, dim(ExprM.RawCounts.filter)[1]-(area1+area2-n12)),  nrow = 2)
  cluster.overlap
  ft.pvalue = fisher.test(cluster.overlap)$p.value
  
  # .3 VennDiagram
  pdf(paste("figures/", exprimentID,"_",rare.cluster,"_diff_gene_overlap_Gini.pdf", sep=""), height = 6, width = 6, useDingbats = FALSE)
  overlapGenes <- intersect(cluster.diffgene, GeneList.final)
  grid.newpage()
  draw.pairwise.venn(length(cluster.diffgene),
                     length(GeneList.final), 
                     length(overlapGenes),
                     category = c(paste("DiffGenes\n", area1), paste("HighGiniGenes\n",area2,"\n","p=",formatC(ft.pvalue, format = "e", digits = 2))), 
                     lty = rep("blank", 2),
                     fill = c("blue","yellow"), 
                     alpha = rep(0.5, 2), 
                     cat.pos = c(0, 0),
                     cat.dist = rep(0.025, 2))
  dev.off()
}

#for common groups, compare to Fano
for(common.cluster in commonGroups){
  # .1 vennplot plot 
  if(file.exists(paste("results/",common.cluster,"_MIST_Sig.RData", sep=""))) {
    load(paste("results/",common.cluster,"_MIST_Sig.RData", sep=""))
    # .2 overlap fisher.test
    rownames(cluster_lrTest.table) = cluster_lrTest.table$Gene
    cluster_lrTest.table$log2fold_change = as.numeric(as.character(cluster_lrTest.table$log2fold_change)) 
    cluster.diffgene <- rownames(cluster_lrTest.table[cluster_lrTest.table$p_value < lr.p_value_cutoff  & abs(cluster_lrTest.table$log2fold_change) > diff.cutoff,])
    length(cluster.diffgene)
    length(intersect(cluster.diffgene, maxFano)) 
    intersect(cluster.diffgene, maxFano)
    
    area1 = length(cluster.diffgene)
    area2 = length(maxFano)
    n12   = length(intersect(cluster.diffgene, maxFano))
    cluster.overlap <- matrix(c(n12, area1-n12 , area2-n12, dim(ExprM.RawCounts.filter)[1]-(area1+area2-n12)),  nrow = 2)
    cluster.overlap
    ft.pvalue = fisher.test(cluster.overlap)$p.value
    
    # .3 VennDiagram
    pdf(paste("figures/", exprimentID,"_",common.cluster,"_diff_gene_overlap_Fano.pdf", sep=""), height = 6, width = 6, useDingbats = FALSE)
    overlapGenes <- intersect(cluster.diffgene, maxFano)
    grid.newpage()
    draw.pairwise.venn(length(cluster.diffgene),
                       length(maxFano), 
                       length(overlapGenes),
                       category = c(paste("DiffGenes\n", area1,"\n","p=",formatC(ft.pvalue, format = "e", digits = 2)), paste("HighFanoGenes\n",area2)), 
                       lty = rep("blank", 2),
                       fill = c("blue","yellow"), 
                       alpha = rep(0.5, 2), 
                       cat.pos = c(0, 0),
                       cat.dist = rep(0.025, 2))
    dev.off()
  }
}