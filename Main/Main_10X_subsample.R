#fixed parameters for GiniClust2
minCellNum           = 3                                                # filtering, remove genes expressed in fewer than minCellNum cells
minGeneNum           = 2000                                             # filtering, remove cells expressed in fewer than minGeneNum genes
expressed_cutoff     = 1                                                # filtering, for raw counts
gini.bi              = 0                                                # fitting, default is 0, for qPCR data, set as 1. 
log2.expr.cutoffl    = 0                                                # cutoff for range of gene expression   
log2.expr.cutoffh    = 20                                               # cutoff for range of gene expression 
Gini.pvalue_cutoff   = 0.0001                                           # fitting, Pvalue, control how many Gini genes chosen
Norm.Gini.cutoff     = 1                                                # fitting, NormGini, control how many Gini genes chosen, 1 means not used.
span                 = 0.9                                              # parameter for LOESS fitting
outlier_remove       = 0.75                                             # parameter for LOESS fitting
GeneList             = 1                                                # parameter for clustering, 1 means using pvalue, 0 means using HighNormGini
Gamma                = 0.9                                              # parameter for clustering
diff.cutoff          = 1                                                # MAST analysis, filter genes that don't have high a log2_foldchange to reduce gene num
lr.p_value_cutoff    = 1e-5                                             # MAST analysis, pvalue cutoff to identify differentially expressed genes
CountsForNormalized  = 100000                                           # if normalizing- by default not used
Rfundir              = "/Users/Daphne/Documents/Yuan/GiniClust2_V1/Rfunction/"     
                                                                        # where GiniClust2 R functions are stored

#dataset specific parameters:
MinPts               = 3                                                # parameter for DBSCAN
eps                  = 0.45                                             # parameter for DBSCAN
mycols               = c("grey50","greenyellow","red","blue","black","green","orange","purple","yellow","navy","magenta")
                                                                        # color setting for tSNE plot
perplexity_G         = 30                                               # parameter for Gini tSNE
perplexity_F         = 30                                               # parameter for Fano tSNE
max_iter_G           = 1000                                             # parameter for Gini tSNE
max_iter_F           = 1000                                             # parameter for Fano tSNE
ks                   = c(2,2,2,3,3,3,3)                                 # a range of k's for k-means for subsampled data depending on rarity: use k=2 for rarer, k=3 for more common
gap_statistic        = FALSE                                            # whether the gap statistic should be used to determine k
K.max                = 10                                               # if using the gap statistic, highest k that should be considered
automatic_eps        = TRUE                                             # whether to determine eps using KNN- for consistency we use the same eps as full data set here
automatic_minpts     = TRUE                                             # whether to determine MinPts based on the size of the data set                                          
workdir              = "/Users/Daphne/Documents/Yuan/GiniClust2_V1/Proj/10X_subsampled/"     
                                                                        # where you put the data and results


setwd(workdir)
dir.create(file.path(workdir, "results"), showWarnings = FALSE) #folder to save results
dir.create(file.path(workdir, "figures"), showWarnings = FALSE) #folder to save figures
#load packages and functions

source(paste(Rfundir,"GiniClust2_packages.R",sep=""))
source(paste(Rfundir,"GiniClust2_functions.R",sep=""))

#generate 140 data sets, with different proportions
source(paste(Rfundir,"Generate_10X_datasets.R",sep=""))

#for each of 140 data sets, run GiniClust2
#for each, plot a barplot comparing the reference and the GiniClust2 result
for (j in 1:7){
  for (i in 1:20){
    k=ks[j]
    exprimentID<-paste("10X_rare",j,"_",i,sep="")
    source(paste(Rfundir,"GiniClust2_fitting.R",sep=""))
    ExprM.RawCounts<-ExprM.RawCounts.filter
    source(paste(Rfundir,"GiniClust2_Gini_clustering.R",sep=""))
    print(table(P_G))
    source(paste(Rfundir,"GiniClust2_Fano_clustering.R",sep=""))
    print(table(P_F))
    source(paste(Rfundir,"GiniClust2_consensus_clustering.R",sep=""))
    print(table(finalCluster))
    
    #plot
    names<-colnames(ExprM.RawCounts.filter)
    barData4<-cbind(finalCluster,names)
    colnames(barData4)<-c("clustering","standard")
    print(ggplot(data=as.data.frame(barData4), aes(x=as.factor(standard),fill=as.factor(clustering),xlab="Developmental Stage")) +
            geom_bar(stat="count")+labs(x="Cell Type",y="Count",fill="Clusters"))
    ggsave(paste("figures/", exprimentID, "GiniClust2_vs_Reference.pdf", sep=""))
  }
}

