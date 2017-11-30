#Preprocess inDrop day 4 post-LIF data

ExprM.RawCounts <- read.delim("data/GSM1599498_ES_d4_LIFminus.csv", sep=",", head=F)

#raw data
title=c("Symbol");
for(i in 2:ncol(ExprM.RawCounts)){
  title=c(title,paste(exprimentID, ".Cell_",i-1,sep=""))
};
colnames(ExprM.RawCounts)=title
rownames(ExprM.RawCounts)=ExprM.RawCounts[,1]
ExprM.RawCounts= ExprM.RawCounts[,-1]
colsum <- apply(ExprM.RawCounts,2,sum)
ExprM.normCounts <- t(t(ExprM.RawCounts)*CountsForNormalized/colsum)
write.table(ExprM.RawCounts, file=paste("results/", exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(ExprM.normCounts, file=paste("results/", exprimentID, "_normCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
save(ExprM.RawCounts, ExprM.normCounts, file=paste("results/",exprimentID, "_ExprM.RData",sep=""))
