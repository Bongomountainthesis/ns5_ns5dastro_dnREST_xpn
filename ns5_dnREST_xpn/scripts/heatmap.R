#!/usr/local/bin/Rscript

options(stringsAsFactors=FALSE)

library(beadarray)
library(gplots)
library(RColorBrewer)

load("results/BSData.quantile.RData")

limma <- read.csv("results/significantuniquegenes.csv")

E <- exprs(BSData.quantile)

limma.o <- limma[order(limma[,"adj.P.Val"],decreasing = FALSE),]
top.ids <- limma.o[1:500,c("ID","EnsemblID")]
#sig.ids <- limma[which(limma[,"adj.P.Val"] <= 0.000001),]

ids <- as.character(top.ids[,"ID"])

filteredE <- E[ids,]

rownames(filteredE) <- top.ids[,"EnsemblID"]

### average filteredE

ns.cols <- grep("NS", colnames(filteredE))
n.cols  <- grep("N[^S]", colnames(filteredE)) 

aves <- apply(filteredE, 1, function(x){  
    NS.ave  <- sum(x[ns.cols])/length(ns.cols)
    N.ave   <- sum(x[n.cols])/length(n.cols)
    return(c(NS.ave, N.ave))
 }
)

aves <- t(aves)
colnames(aves)<-c("NS.ave", "N.ave")

##remove NAs

aves <- aves[which(!(is.na(rownames(aves))),]

distance <- dist(aves)
cluster <- hclust(distance, method = "ward")
dendrogram <- as.dendrogram(cluster)

## then get order from dendrogram

dendro_order <- rev(order.dendrogram(dendrogram))

aves_dendro <- aves[dendro_order,]

col.pal <- rev(brewer.pal(11,"RdYlBu")[c(1,2,3,4,5,7,8,9,10,11)])

postscript(file="results/heatmap_avg.ps", horizontal=FALSE)
heatmap.2(aves_dendro,
		Colv=NA,
		Rowv=NA,
		labRow="",
		col=col.pal,
		scale="none",
		key=TRUE,
		keysize=0.75,
		symkey=FALSE,
		density.info="none",
		trace="none", 
		labCol=c("NSC","Neuron"),
		colsep=1:2,
		sepcolor="white",
		sepwidth=0.1
	)
dev.off()

###now take blocks of colour and assign to GO terms for figure

#split here
#ENSMUSG00000024109
#ENSMUSG00000048349
#ENSMUSG00000061758

#add column with numbers

aves_dendro <- cbind(aves_dendro,seq(1,nrow(aves_dendro),1))

go1 <- rownames(aves_dendro[1:91,])
go2 <- rownames(aves_dendro[92:180,])
go3 <- rownames(aves_dendro[182:301,])
go4 <- rownames(aves_dendro[302:495,])

library(gtools)
go.ids <- smartbind(go1,go2,go3,go4)
go.ids <- t(go.ids)

write.csv(go.ids, file = "results/go_ids_for_heatmap_clusters.csv")





















