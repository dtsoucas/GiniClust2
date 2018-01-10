#Find differentially expressed genes for each cluster in "finalCluster"

load(paste("results/", exprimentID,"_ExprM.filter.RData", sep=""))  #"ExprM.RawCounts.filter"  "ExprM.normCounts.filter"
load(file=paste("results/", exprimentID,"_FinalClustering.RData",sep=""))

pseudo.count = 0.1
data.used.log2   <- log2(ExprM.RawCounts.filter+pseudo.count)
data.used.log2<-as.matrix(data.used.log2)
colnames(data.used.log2)<-c(1:dim(data.used.log2)[2])

for(current.cluster in unique(finalCluster)){
  
  cells.symbol.list2     = colnames(data.used.log2)[which(finalCluster==current.cluster)]
  cells.coord.list2      = match(cells.symbol.list2, colnames(data.used.log2))                          
  cells.symbol.list1     = colnames(data.used.log2)[which(finalCluster != current.cluster)]
  cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.log2))   
  data.used.log2.ordered  = cbind(data.used.log2[,cells.coord.list1], data.used.log2[,cells.coord.list2])
  group.v <- c(rep(0,length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
  #ouput
  log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
  Auc <- m.auc(data.used.log2.ordered, group.v)
  bigtable <- data.frame(cbind(log2.stat.result, Auc))
  
  diff.cutoff = 1
  DE <- bigtable[bigtable$log2_fc >diff.cutoff,] 
  dim(DE)
  if (dim(DE)[1]>0) {
    data.1                 = data.used.log2[,cells.coord.list1]
    data.2                 = data.used.log2[,cells.coord.list2]
    genes.list = rownames(DE)
    log2fold_change        = cbind(genes.list, DE$log2_fc)
    colnames(log2fold_change) = c("gene.name", "log2fold_change")
    counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
    groups  = c(rep("Cluster_Other", length(cells.coord.list1) ), rep(current.cluster, length(cells.coord.list2) ) )
    groups  = as.character(groups)
    data_for_MIST <- as.data.frame(cbind(rep(rownames(counts), dim(counts)[2]), melt(counts),rep(groups, each = dim(counts)[1]), rep(1, dim(counts)[1] * dim(counts)[2]) ))
    colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
    vbeta = data_for_MIST
    vbeta.fa <-FromFlatDF(vbeta, idvars=c("Subject.ID"),
                          primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                          geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                          phenovars=c('Population'), id='vbeta all')
    vbeta.1 <- subset(vbeta.fa,Number.of.Cells==1)
    # .3 MAST 
    head(colData(vbeta.1))
    zlm.output <- zlm(~ Population, vbeta.1, method='bayesglm', ebayes=TRUE)
    show(zlm.output)
    coefAndCI <- summary(zlm.output, logFC=TRUE)
    zlm.lr <- lrTest(zlm.output, 'Population')
    zlm.lr_pvalue <- melt(zlm.lr[,,'Pr(>Chisq)'])
    zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]
    
    lrTest.table <-  merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
    colnames(lrTest.table) <- c("Gene", "test.type", "p_value", paste("log2.mean.", "Cluster_Other", sep=""), paste("log2.mean.",current.cluster,sep=""), "log2fold_change", "Auc")
    cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)),]
    
    #. 4 save results
    write.csv(cluster_lrTest.table, file=paste("results/",current.cluster,"_lrTest_Sig.csv", sep=""))
    save(cluster_lrTest.table, file=paste("results/",current.cluster,"_MIST_Sig.RData", sep=""))
  }
  
}
