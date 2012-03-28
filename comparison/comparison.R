

options(stringsAsFactors = FALSE)

ns5 <- read.csv(file="ns5_dnREST_xpn/results/ns5_limma_results.csv")

astro <- read.csv(file="ns5dastro_dnREST_xpn/results/astro_limma_rd.csv")

##tidy up astro file

#astro <- astro[,c(5,25,22,1,2,3,4,6,7,8,9,10,11)]

######### astros - need to tidy up astro table as its a state...

astro.symbol <- as.character(astro[,"Proportion_UCSC_transcripts"])

astro.ensembl <- as.character(astro[,"Proportion_Ensembl_transcripts"])

astro.symbol.brackets <- substr(gsub('.*\\(','',astro.symbol),0,nchar(gsub('.*\\(','',astro.symbol))-1)

astro.ensembl.brackets <- substr(gsub('.*\\(','',astro.ensembl),0,nchar(gsub('.*\\(','',astro.ensembl))-1)

#make astro dataframe and save it this time...

astro.tidy <- astro[,c(1,2,3,4,6,9,10)]

astro.res <- cbind(astro.tidy,astro.symbol.brackets,astro.ensembl.brackets)

astro.retidy <- astro.res[,c(9,8,1,2,3,4,5,6,7)]

colnames(astro.retidy) <- c("EnsemblID","Symbol","Chr","Start","End","Strand","logFC","P.Value","adj.P.Val")

write.csv(astro.retidy, file = "ns5_dnREST_xpn/results/ns5dastro_dnREST_xpn_tidy.csv")

## remove duplicates

astro.o <- astro.retidy[order(abs(astro.retidy[,"logFC"]),decreasing = TRUE),]
astro <- astro.o[!duplicated(astro.o[,"EnsemblID"],decreasing = TRUE),]

ns5 <- ns5[order(ns5[,"adj.P.Val"], decreasing = FALSE),]
ns5 <- ns5[!duplicated(ns5[,"EnsemblID"]),]

##remove NAs

astro <- astro[which(!(is.na(astro[,"EnsemblID"]))),]
ns5 <- ns5[which(!(is.na(ns5[,"EnsemblID"]))),]

##remove genes from astro set without EnsemblID - blank space

astro <- astro[grep("\\w",astro[,"EnsemblID"]),]

########### now can use the data

##take all genes at adj.p.val of 0.05

ns5 <- ns5[which(ns5[,"adj.P.Val"] <= 0.05),]
astro <- astro[which(astro[,"adj.P.Val"] <= 0.05),]

##up and down regulated
ns5.up <- ns5[which(ns5[,"logFC"]>=1),"EnsemblID"]
ns5.down <- ns5[which(ns5[,"logFC"]<=-1),"EnsemblID"]

astro.up <- astro[which(astro[,"logFC"]>=1),"EnsemblID"]
astro.down <- astro[which(astro[,"logFC"]<=-1),"EnsemblID"]

##write out for GO terms

write.csv(ns5, file = "ns5_dnREST_sig_changing.csv")
write.csv(astro, file = "astro_dnREST_sig_changing.csv")

#### now work out shared/unique genes

ns5.ids <- ns5[,"EnsemblID"]

astro.ids <- astro[,"EnsemblID"]

shared.ids <- intersect(ns5.ids,astro.ids)

#### 40 genes shared

##stick them into dataframe....

ns5.shared <- ns5[which(ns5[,"EnsemblID"] %in% shared.ids),]

astro.shared <- astro[which(astro[,"EnsemblID"] %in% shared.ids),]

shared.merge <- merge(ns5.shared,astro.shared, by.x = "EnsemblID", by.y = "EnsemblID", suffixes = c("_NS5","_Astro"))

##plot changes in FC

postscript(file = "results/shared_gene_expression_changes_after_dnREST", horizontal = FALSE)
plot(shared.merge[,"logFC_NS5"],shared.merge[,"logFC_Astro"],
		pch = 20,
		main = "Shared Gene Expression Changes after dnREST",
		xlab = "NS5 Gene Expression Changes",
		ylab = "Astrocyte Gene Expression Changes"
		)
dev.off()


#tidy and save

shared.merge.tidy <- shared.merge[,c(1,2,3,5,6,7,8,9,10,13,14,21,22,23)]

write.csv(shared.merge.tidy, file = "gene_changes_shared_between_ns5_and_astro_treated_with_DNREST.csv")

###### find unique genes now...

ns5.unique <- ns5.ids[which(!(ns5.ids %in% astro.ids))]

astro.unique <- astro.ids[which(!(astro.ids %in% ns5.ids))]

# pull out data from dataframes

ns5.unique.df <- ns5[which(ns5[,"EnsemblID"] %in% ns5.unique),]

astro.unique.df <- astro[which(astro[,"EnsemblID"] %in% astro.unique),]

write.csv(ns5.unique.df, file = "gene_changes_unique_to_ns5_treated_with_DNREST.csv")

write.csv(astro.unique.df, file = "gene_changes_unique_to_astro_treated_with_DNREST.csv")

## up shared

up.shared <- intersect(ns5.up, astro.up)
down.shared <- intersect(ns5.down, astro.down)

## no shared down genes - get shared up

up.shared.df <- shared.merge[which(shared.merge[,"EnsemblID"] %in% up.shared),"symbol"]












