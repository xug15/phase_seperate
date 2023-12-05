
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")

library(edgeR)
#set the work directory.
setwd('D:/a-document/github/phase_seperate/data/a9-genecount')
#read the gene counts.
genecount=read.table('merge.gene.tsv',header = T,sep = "\t",row.names = c(1))
head(genecount)
#mutation vs wild type.
mut.wt=genecount[,c(1,2,3,7,8,9)]
#over expression vs wild type.
oe.wt =genecount[,c(4,5,6,7,8,9)]
#normolize to one million counts.
mut.wt.cpm=cpm(mut.wt)
oe.wt.cpm=cpm(oe.wt)





