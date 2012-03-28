#!/usr/local/bin/Rscript

library(beadarray)
library(limma)
library(biomaRt)
library(illuminaMousev2BeadID.db)

		#load the data
BSData <- get(load("data/BSData.quantile.RData"))
E <- exprs(BSData)

design<-matrix(0,nrow=(ncol(E)), ncol=3)
colnames(design) <- c("E11.5","E14.5","P1")
rownames(design) <- colnames(E)
design[1:2,1] <- 1
design[3:6,2] <- 1
design[7:8,3] <- 1
cont.matrix<-makeContrasts(E11.5vE14.5=E14.5-E11.5, E14.5vP1=P1-E14.5, levels=design)


fit<-lmFit(E, design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

ids = rownames(E)
#remove ILMN_ prefix
ids = substr(ids, 6,nchar(ids))

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

write.fit(ebFit, file="results/limma_ebfit11.5_14.5.csv", adjust="BH")
data<-read.table("results/limma_ebfit11.5_14.5.csv", sep="\t", header=T)


new.data<- topTable(ebFit, number=nrow(E))
rownames(new.data)<-new.data$ID
new.data<-new.data[order(new.data[,"P.Value"]),]
write.csv(new.data,"results/limma_results11.5_14.5.csv",row.names=F)


