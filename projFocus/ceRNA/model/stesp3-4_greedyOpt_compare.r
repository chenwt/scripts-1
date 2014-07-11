rm(list=ls())

filem = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul8wCNV/optCorr.result_flex.test.summary.tsv"
filewCnv = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul8wCNV/optCorr.result_flex_maxWcnv.test.summary.tsv"

data1 = read.table(filem, sep="\t", header = T,stringsAsFactors=F)
data2 = read.table(filewCnv, sep = "\t" , header = T, stringsAsFactors=F)

ji = sapply(1:nrow(data1), FUN=function(x) {
  mut = unlist(strsplit(data1$optGene.optSmp[x], ";"))
  cnv = unlist(strsplit(data2$optGene.optSmp[x], ";"))
  cntm = length(mut); cntc = length(cnv)
  cntmatch = length(intersect(mut,cnv))
  cntmonly = length(setdiff(mut,cnv))
  cntcnvonly = length(setdiff(cnv, mut))
  cntunion = length(unique(c(mut,cnv)))
  return(c(cntm, cntc, cntmatch, cntmonly, cntcnvonly, cntunion, cntmatch/ cntm * 100, cntmatch/cntunion * 100))
} )

ji = t(ji)
ji
par(mfrow=c(2,2))
hist(ji[,7], col = "lightblue")
boxplot(ji[,1], col = 'lightblue')
boxplot(ji[,2], col = "lightblue")
boxplot(ji[,3], col="lightblue")

plot(sort(ji[,1]))
plot(density(ji[,3]))

t.test(ji[,2], ji[,1],alternative='greater')
