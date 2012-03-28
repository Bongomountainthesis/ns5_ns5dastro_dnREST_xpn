#!/usr/local/bin/Rscript

library(beadarray)


dataFile = "data/Raw_Data_UnNorm_Rewritten.txt"



BSData <- readBeadSummaryData(dataFile=dataFile, 
                              skip=0, 
				ProbeID="ProbeID",
                              columns = list(exprs = "AVG_Signal", 
                                se.exprs="BEAD_STDERR", 
                                NoBeads = "Avg_NBEADS", 
                                Detection="Detection Pval"
					     ),
                              )

labels <- read.csv("data/labels.tsv", sep="\t", as.is=T)
rownames(labels) <- labels$Samples

#BSData<-BSData[,1:8]
BSData<-BSData[,c(12,4,23,11,9,15,22,7)]
labels <- labels[1:8,]
sampleNames(BSData)<-labels[sampleNames(BSData),"Label"]

#fix.col<-function(x){
#	x[x<=0]<-1 
#	return(x)
#}

#exprs(BSData) <- apply(exprs(BSData),2,fix.col) 

save(BSData, file="data/BSData.RData")

postscript(file="results/Boxplotprenorm.ps", horizontal=FALSE)
boxplot(as.data.frame(log2(exprs(BSData))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,4,4,4,4,6,6), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYprenorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData), arrays = 1:8,  pch = 16)
dev.off()

postscript(file="results/plotDENSITYprenorm.ps", horizontal=FALSE)
E <- log2(exprs(BSData))
plot(density(E[,1]))
for(i in 1:8){
  lines(density(E[,i]),col=i)
}
dev.off()



BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")


postscript(file="results/Boxplotpostnorm.ps", horizontal=FALSE)
boxplot(as.data.frame((exprs(BSData.quantile))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,4,4,4,4,6,6), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYpostnorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData.quantile), arrays = 1:8, log = FALSE, pch = 16)
dev.off()

save(BSData.quantile, file="data/BSData.quantile.RData")

postscript(file="results/plotDENSITYpostnorm.ps", horizontal=FALSE)
E <- exprs(BSData.quantile)
plot(density(E[,1]))
for(i in 1:8){
  lines(density(E[,i]),col=i)
}
dev.off()


