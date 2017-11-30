#Preprocess simulated data

ExprM.RawCounts <- read.table("data/Data_2000_1000_10_6_4_3.xls", sep="\t", head=TRUE, row.names=1) 
dim(ExprM.RawCounts)

colsum <- apply(ExprM.RawCounts,2,sum)
ExprM.normCounts <- t(t(ExprM.RawCounts)*CountsForNormalized/colsum)
write.table(ExprM.RawCounts, file=paste("results/", exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(ExprM.normCounts, file=paste("results/", exprimentID, "_normCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
save(ExprM.RawCounts, ExprM.normCounts, file=paste("results/",exprimentID, "_ExprM.RData",sep=""))
dim(ExprM.RawCounts)
