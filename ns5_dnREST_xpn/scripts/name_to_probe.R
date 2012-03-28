#!/usr/local/bin/Rscript

require(illuminaMousev2BeadID.db)
require(illuminaMousev2.db)


args <- commandArgs(trailingOnly=TRUE)
symbol <- args[1]

annot <- mget(ls(illuminaMousev2SYMBOL), illuminaMousev2SYMBOL, ifnotfound = NA)

res <- annot[grep(symbol, annot, ignore.case=T)]

if(length(res)>0){
res
}else{stop("No matches found")}
