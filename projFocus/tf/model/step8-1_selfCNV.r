rm(list=ls())
setwd("/Users/jh3283/projFocus/07252014/step1_selfCNV/")

cnvfile = "FOXC1.cnv"
expfile = "FOXC1.exp"

dcnv = read.table(cnvfile,header=T,stringsAsFactors=F)
dexp = read.table(expfile, header=T, stringsAsFactors=F)
tgene = rownames(dexp)

###save this part later---

require(useful)
commSmps = intersect(colnames(dcnv), colnames(dexp))

dcnv = subset(dcnv, select=commSmps) ; dexp = subset(dexp, select=commSmps)
rownames(dcnv) <- "selfcnv"
rownames(dexp) <- "selfexp"


###prepare data
dmodel = data.frame(t(rbind(dexp, dcnv)))
fitlm = lm(selfexp ~ selfcnv, dmodel)

save(tgene, fitlm, file="1_selfcnv.rda")
