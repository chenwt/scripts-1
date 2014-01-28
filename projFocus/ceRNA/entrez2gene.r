#!/usr/bin/Rscript
#use lab data: /ifs/scratch/c2b2/ac_lab/jh3283/ref/aclab/entrez2gene.rda

setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/knowledgeBase")
load("/ifs/scratch/c2b2/ac_lab/jh3283/ref/aclab/entrez2gene.rda")
input="/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/knowledgeBase/GWAS_catalog_brca_allGene.txt"
id = vapply(unlist(read.table(input)),FUN=as.character,'a')
name = sort(unique(entrez2gene[id]))
name = as.matrix(na.omit(name))
write.table(name,"/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/knowledgeBase/GWAS_catalog_brca_allGeneName.txt",quote=F,row.names=F,col.names=F)


