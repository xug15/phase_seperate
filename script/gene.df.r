
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

group <- factor(c(2,2,2,1,1,1))

y <- DGEList(counts=mut.wt,group=group)

keep <- filterByExpr(y)

y <- y[keep,,keep.lib.sizes=FALSE]

y <- calcNormFactors(y)

design <- model.matrix(~group)

y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)


fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
result=topTags(lrt,n=1000000)

result
#lrt$table

lrt$fitted.values

resultt=result$table
#head(resultt)
up=resultt[resultt$logFC>1 & resultt$PValue<0.05,]
do=resultt[resultt$logFC< -1 & resultt$PValue<0.05,]
head(up)
head(do)
write.table(up,"mut.vs.wt.up.tsv",quote=F,sep="\t")
write.table(do,"mut.vs.wt.down.tsv",quote=F,sep="\t")
write.table(mut.wt.cpm,"mut.vs.wt.exp.tsv",quote=F,sep="\t")
###

y <- DGEList(counts=oe.wt,group=group)

keep <- filterByExpr(y)

y <- y[keep,,keep.lib.sizes=FALSE]

y <- calcNormFactors(y)

design <- model.matrix(~group)

y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)


fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
result=topTags(lrt,n=1000000)

result
#lrt$table

lrt$fitted.values

resultt=result$table
#head(resultt)
up=resultt[resultt$logFC>1 & resultt$PValue<0.05,]
do=resultt[resultt$logFC< -1 & resultt$PValue<0.05,]
head(up)
head(do)
write.table(up,"oe.vs.wt.up.tsv",quote=F,sep="\t")
write.table(do,"oe.vs.wt.down.tsv",quote=F,sep="\t")
write.table(mut.wt.cpm,"oe.vs.wt.exp.tsv",quote=F,sep="\t")


