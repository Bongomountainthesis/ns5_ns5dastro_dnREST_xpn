options(stringsAsFactors = FALSE)

library(beadarray)
library(limma)
library(biomaRt)

library(illuminaMousev2BeadID.db) 

library(illuminaMousev2.db) 


dataFile = "data/Unnormaliseddata_MBurney_280212.txt"

BSData <- readBeadSummaryData(dataFile=dataFile,
				skip=0,
				ProbeID="ProbeID",
				sep = "\t"
				)

save(BSData, file = "data/BSData.RData")

##QC plots

postscript(file="results/Boxplotprenorm.ps", horizontal=FALSE)
boxplot(as.data.frame(log2(exprs(BSData))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,4,2,4,2,4), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYprenorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData), arrays = 1:6,  pch = 16)
dev.off()

postscript(file="results/plotDENSITYprenorm.ps", horizontal=FALSE)
E <- log2(exprs(BSData))
plot(density(E[,1]))
for(i in 1:6){
  lines(density(E[,i]),col=i)
}
dev.off()

# Quantile Normalise
BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")


postscript(file="results/Boxplotpostnorm.ps", horizontal=FALSE)
boxplot(as.data.frame((exprs(BSData.quantile))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,4,2,4,2,4), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYpostnorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData.quantile), arrays = 1:6, log = FALSE, pch = 16)
dev.off()

save(BSData.quantile, file="data/BSData.quantile.RData")

postscript(file="results/plotDENSITYpostnorm.ps", horizontal=FALSE)
E <- exprs(BSData.quantile)
plot(density(E[,1]))
for(i in 1:6){
  lines(density(E[,i]),col=i)
}
dev.off()

## LIMMA analysis

##rearrange E

E <- exprs(BSData.quantile)

E <- E[,c(3,1,2,6,5,4)]

### NS5 treated with empty vector and with DomNeg REST

design<-matrix(0,nrow=(ncol(E)), ncol=2)

colnames(design) <- c("Vector","DomNeg")
rownames(design) <- colnames(E)

design[c(1,2,3),1] <- 1
design[c(4,5,6),2] <- 1

cont.matrix <- makeContrasts(VectorvsDomNeg = DomNeg-Vector,
				levels = design
				)


fit<-lmFit(E, design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

ids = rownames(E)

symbol <- mget(ids, illuminaMousev2BeadIDSYMBOL, ifnotfound = NA)
sum(sapply(symbol, length) > 1)
symbol <- as.character(symbol)
length(which(symbol=="NA"))

ensembl = mget(ids, illuminaMousev2BeadIDENSEMBL, ifnotfound = NA)
length(which(is.na(ensembl)))

sum(sapply(ensembl, length) > 1)
crosshyb <- which(( sapply(ensembl, length) ) > 1)
length(crosshyb)
ensembl[crosshyb] <- NA
ensembl <- as.character(ensembl)
ensembl[ensembl=="NA"] = NA
length(which(is.na(ensembl)))

ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")


filters <- "ensembl_gene_id"
values <- ensembl[!is.na(ensembl)]
attributes <- c("ensembl_gene_id", "chromosome_name","start_position", "end_position", "strand", "description")
ens.anno <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)
rownames(ens.anno)<-ens.anno[,1]

anno <- data.frame(
                   ID = as.character(ids),
                   EnsemblID = ensembl,
                   symbol=symbol,
                   ens.anno[ensembl,],
                   stringsAsFactors=F
              )
rownames(anno) <- anno[,"ID"]

ebFit$genes = anno

write.fit(ebFit, file="results/limma_ebfit.csv", adjust="BH")
data<-read.table("results/limma_ebfit.csv", sep="\t", header=T)

new.data<- topTable(ebFit, number=nrow(E))
rownames(new.data)<-new.data$ID
new.data<-new.data[order(new.data[,"adj.P.Val"],decreasing = FALSE),]
write.csv(new.data,"results/limma_results.csv",row.names=F)

###plot some heatmap for Noel

#take most significantly changing genes and plot on heatmap

top <- new.data[!duplicated(new.data[,"EnsemblID"]),]

top <- top[which(!(is.na(top[,"EnsemblID"]))),]

top <- top[1:1000,]

##cut off on pval
#top <- top[which(top[,"adj.P.Val"] <= 0.05),]

top.probes <- top[,"ID"]

topE <- E[top.probes,]

##average probes for each sample

aves <- apply(topE, 1, function(x){  
    NS.vec  <- sum(x[1:3])/3
    NS.dneg   <- sum(x[4:6])/3
    return(c(NS.vec, NS.dneg))
 }
)

aves <- t(aves)
colnames(aves)<-c("NS.vec", "N.dneg")

distance <- dist(aves)
cluster <- hclust(distance, method = "ward")
dendrogram <- as.dendrogram(cluster)

dendro_order <- rev(order.dendrogram(dendrogram))

aves_dendro <- aves[dendro_order,]

library(gplots)
library(RColorBrewer)

postscript(file = "results/heatmap.eps", horizontal = FALSE)
heatmap.2(aves,
	Colv = NA, 
	labRow = NA,
	scale = "row",
	density.info = "none",
	trace = "none",
	col=brewer.pal(9, "Blues"),
	colsep = 1:2,
	sepwidth = 0.05,
	sepcol = "white",
	keysize = 0.75
	)

dev.off()

