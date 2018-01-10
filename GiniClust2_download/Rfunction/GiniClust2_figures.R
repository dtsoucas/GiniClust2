#composite tSNE visualization

#load Gini and Fano dimensions
Rtsne_map_G<-read.table(paste("results/", exprimentID,"_Rtnse_Gini_coord2.csv", sep=""), sep=",")
Rtsne_map_F<-read.table(paste("results/", exprimentID,"_Rtnse_Fano_coord2.csv", sep=""), sep=",")

n<-length(unique(finalCluster))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


#load final clustering ids
load(file=paste("results/", exprimentID,"_FinalClustering.RData",sep=""))

#load filtered data set
load(paste("results/", exprimentID,"_ExprM.filter.RData", sep=""))  #"ExprM.RawCounts.filter"  "ExprM.normCounts.filter"


#plot visualizations

mydata<-cbind(Rtsne_map_F[,1],Rtsne_map_F[,2],Rtsne_map_G[,1])
colnames(mydata)<-c("1m","2m","3m")
mydata<-as.data.frame(mydata)
pdf(paste("figures/Figures_", exprimentID, "_tsne_plot_combined.pdf", sep=""))
with(mydata, {
   scatterplot3d(Rtsne_map_F[,1],   # x axis
                 Rtsne_map_F[,2],     # y axis
                 Rtsne_map_G[,1], color=mycols[finalCluster], pch=19,cex.symbols=0.4,xlab="Fano 1",ylab="Fano 2",zlab="Gini",cex.axis=1.2)
})
legend("topright",col=mycols[1:length(unique(finalCluster))],legend=(1:length(unique(finalCluster))),lty=c(1,1), lwd=2,bg="white",cex=0.5)
dev.off()


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
