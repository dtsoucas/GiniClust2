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
diff.cutoff          = 1                                                # MAST analysis, filter genes that don't have high log2_foldchange to reduce gene num
lr.p_value_cutoff    = 1e-5                                             # MAST analysis, pvalue cutoff to identify differentially expressed genes
CountsForNormalized  = 100000						# if normalizing- by default not used
Rfundir              = "/Users/Daphne/Documents/Yuan/GiniClust2_V1/Rfunction/"
									# where GiniClust2 R functions are stored

#dataset-specific parameters:
MinPts               = 3                                                # parameter for DBSCAN
eps                  = 0.34                                             # parameter for DBSCAN
mycols               = c("grey50","greenyellow","red","blue","black","orange")                                 
									# color setting for tSNE plot
perplexity_G         = 30						# parameter for Gini tSNE
perplexity_F         = 30						# parameter for Fano tSNE
max_iter_G           = 1000 						# parameter for Gini tSNE
max_iter_F           = 1000   						# parameter for Fano tSNE
k                    = 2 						# k for k-means step
gap_statistic        = FALSE 						# whether the gap statistic should be used to determine k- here will also yield 2
K.max                = 10 						# if using the gap statistic, highest k that should be considered
automatic_eps        = FALSE 						# whether to determine eps using KNN
automatic_minpts     = FALSE 						# whether to determine MinPts based on the size of the data set
workdir              = "/Users/Daphne/Documents/Yuan/GiniClust2_V1/Proj/inDrop_day4/"                              
									# where you put the data and results
exprimentID          = "d4"                              		# experiment or data set ID

setwd(workdir)
dir.create(file.path(workdir, "results"), showWarnings = FALSE) #folder to save results
dir.create(file.path(workdir, "figures"), showWarnings = FALSE) #folder to save figures
#load packages and functions
source(paste(Rfundir,"GiniClust2_packages.R",sep=""))
source(paste(Rfundir,"GiniClust2_functions.R",sep=""))

#Preprocessing the data
source(paste(Rfundir,"PreProcess_for_inDrop_day4.R",sep=""))
source(paste(Rfundir,"GiniClust2_filtering_RawCounts.R",sep=""))

#Gini-based clustering steps
source(paste(Rfundir,"GiniClust2_fitting.R",sep=""))
source(paste(Rfundir,"GiniClust2_Gini_clustering.R",sep=""))

table(P_G) #P_G is the Gini-based clustering result

source(paste(Rfundir,"GiniClust2_Gini_tSNE.R",sep="")) #visualization of GiniClust results using tSNE

#Fano-based clustering steps
source(paste(Rfundir,"GiniClust2_Fano_clustering.R",sep=""))

table(P_F) #P_F is the Fano-based clustering result

source(paste(Rfundir,"GiniClust2_Fano_tSNE.R",sep="")) #visualization of k-means results using tSNE

#weighted consensus clustering
source(paste(Rfundir,"GiniClust2_consensus_clustering.R",sep=""))

table(finalCluster) #finalCluster is the weighted consensus clustering result

#final analyses
source(paste(Rfundir,"GiniClust2_DE.R",sep="")) #find differentially expressed genes for each finalCluster
source(paste(Rfundir,"GiniClust2_figures.R",sep="")) #plot composite tSNE and gene-overlap venn diagrams

