source('gsea_will.R')
mexp = read.table('C:/Documents and Settings/jh3283/My Documents/Dropbox/1Research/AML/data/target_AML_uniq_probes_177.exp',sep='\t', comment.char='%', as.is=TRUE, header=TRUE)
mexp.matrix = as.matrix(mexp[,3:ncol(mexp)])
rownames(mexp.matrix) = mexp[,2]

phyne <- read.table("C:/Documents and Settings/jh3283/My Documents/Dropbox/1Research/AML/data/group.cls",skip = 2, as.is=TRUE,sep="\t")
class(phyne)
index_1 <- unlist(which(phyne==1))-1
index_2 <- unlist(which(phyne==2))-1
index_3 <- unlist(which(phyne==3))-1


nulldist = getNull(mexp=mexp.matrix, pheno1.col=colnames(mexp.matrix)[index_1],pheno2.col=colnames(mexp.matrix)[index_2], np=1000)
geneset = read.table('C:/Documents and Settings/jh3283/My Documents/Dropbox/1Research/AML/data/stirward_hugo.gmt', sep='\t', as.is=TRUE)
geneset.new = unlist(geneset)[2:length(geneset)]
gsea_result = gsea(mexp=mexp.matrix, pheno1.col=colnames(mexp.matrix)[index_1],pheno2.col=colnames(mexp.matrix)[index_2],gs=geneset.new,nullDist=nulldist,plot=T,plotFile="g1_g2.pdf")


nulldist = getNull(mexp=mexp.matrix, pheno1.col=colnames(mexp.matrix)[index_1],pheno2.col=colnames(mexp.matrix)[index_3], np=1000)
geneset = read.table('C:/Documents and Settings/jh3283/My Documents/Dropbox/1Research/AML/data/stirward_hugo.gmt', sep='\t', as.is=TRUE)
geneset.new = unlist(geneset)[2:length(geneset)]
gsea_result = gsea(mexp=mexp.matrix, pheno1.col=colnames(mexp.matrix)[index_1],pheno2.col=colnames(mexp.matrix)[index_3],gs=geneset.new,nullDist=nulldist,plot=T,plotFile="g1_g3.pdf")


nulldist = getNull(mexp=mexp.matrix, pheno1.col=colnames(mexp.matrix)[index_2],pheno2.col=colnames(mexp.matrix)[index_3], np=1000)
gsea_result = gsea(mexp=mexp.matrix, pheno1.col=colnames(mexp.matrix)[index_2],pheno2.col=colnames(mexp.matrix)[index_3],gs=geneset.new,nullDist=nulldist,plot=T,plotFile="g2_g3.pdf")


nulldist = getNull(mexp=mexp.matrix, pheno1.col=colnames(mexp.matrix)[index_1],pheno2.col=colnames(mexp.matrix)[index_2], np=1000)
geneset = read.table('C:/Documents and Settings/jh3283/My Documents/Dropbox/1Research/AML/data/set1allHugo.gmt', sep='\t', as.is=TRUE)
geneset.new = unlist(geneset)[2:length(geneset)]
gsea_result = gsea(mexp=mexp.matrix, pheno1.col=colnames(mexp.matrix)[index_1],pheno2.col=colnames(mexp.matrix)[index_2],gs=geneset.new,nullDist=nulldist,plot=T,plotFile="test.pdf")
