rm(list=ls())
setwd("/Users/jh3283/projFocus/07252014/step1_selfCNV/")

cnvfile = "FOXC1.cdriCNV"
expfile = "FOXC1.exp"

###load data
dexp = read.table(expfile, header = T,stringsAsFactors=F)
dcnv = read.table(cnvfile, header = T, stringsAsFactors=F,sep="\t")

tgene = rownames(dexp); rgene = rownames(dcnv)
###
commSmps = intersect(colnames(dcnv), colnames(dexp))
dcnv = subset(dcnv, select=commSmps) ; dexp = subset(dexp, select=commSmps)
###

dmodel = data.frame(t(rbind(dexp, dcnv)))
colnames(dmodel) <- c(tgene, rgene)

require(glmnet)
fitcv = cv.glmnet(t(dcnv), y = t(dexp),nfolds=10,keep=TRUE)
fitCdriCNV = glmnet(x=t(dcnv), y = t(dexp),lambda=fitcv$lambda.min)

save(tgene,rgene, fitCdriCNV, file="2_cdricnv.rda")
