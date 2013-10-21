myeth <- methyg
dset <- exp
load("../PANCANCER/brca/brca-cnv.rda")

fmyeth <- function(dset, myeth, threshold=1, nn=10, method=c("utest", "ttest", "rea")[3], pval=.05, padj="fdr") {

  #intersection of dset genes and myeth gene
  myeth <- myeth[rownames(myeth) %in% rownames(dset), ]

  myeth <- (abs(myeth)>threshold)*sign(myeth)

  myeth <- filtro.row.matrix(myeth, rowSums(abs(myeth), na.rm=T)>=nn)
  genes <- intersect(rownames(myeth), rownames(dset))
  myeth <- myeth[match(genes, rownames(myeth)), ]
  dset <- dset[match(genes, rownames(dset)), ]
  samples <- intersect(colnames(dset), colnames(myeth))
  myeth <- myeth[, match(samples, colnames(myeth))]
  dset <- dset[, match(samples, colnames(dset))]
  switch(match.arg(method, c("utest", "test", "rea")),
         utest={fmyeth <- myethUtest(dset, myeth)},
         ttest={fmyeth <- myethTtest(dset, myeth)},
         rea={fmyeth <- myethRea(dset, myeth)})
  fmyeth <- p.adjust(pnorm(fmyeth, lower.tail=F), padj)
  return(apply(myeth[fmyeth<pval, ], 1, function(x) x[x!=0]))
}

myethUtest <- function(dset, myeth) {
  myeth[is.na(myeth)] <- 0
  pb <- txtProgressBar(max=nrow(dset), style=3)
  tmp <- sapply(1:nrow(dset), function(i, dset, myeth, pb) {
    zup <- NA
    if (length(which(myeth[i, ]>0))>1) {
      tmp <- wilcox.test(dset[i, myeth[i, ]>0], dset[i, myeth[i, ]<=0])
      zup <- qnorm(tmp$p.value/2, lower.tail=F)*sign(tmp$statistic)
    }
    zdn <- NA
    if (length(which(myeth[i, ]<0))>1) {
      tmp <- wilcox.test(dset[i, myeth[i, ]<0], dset[i, myeth[i, ]>=0])
      zdn <- -qnorm(tmp$p.value/2, lower.tail=F)*sign(tmp$statistic)
    }
    setTxtProgressBar(pb, i)
    return(sum(zup, zdn, na.rm=T)/sqrt(length(which(!is.na(c(zup, zdn))))))
  }, dset=dset, myeth=myeth, pb=pb)
  tmp[is.na(tmp)] <- 0
  names(tmp) <- rownames(dset)
  return(tmp)
}



myethTtest <- function(dset, myeth) {
  myeth[is.na(myeth)] <- 0
  pb <- txtProgressBar(max=nrow(dset), style=3)
  tmp <- sapply(1:nrow(dset), function(i, dset, myeth, pb) {
    zup <- NA
    if (length(which(myeth[i, ]>0))>1) {
      tmp <- t.test(dset[i, myeth[i, ]>0.2], dset[i, myeth[i, ]<=0.2])
      zup <- qnorm(tmp$p.value/2, lower.tail=F)*sign(tmp$statistic)
    }
    zdn <- NA
    if (length(which(myeth[i, ]<0))>1) {
      tmp <- t.test(dset[i, myeth[i, ]<0], dset[i, myeth[i, ]>=0])
      zdn <- -qnorm(tmp$p.value/2, lower.tail=F)*sign(tmp$statistic)
    }
    setTxtProgressBar(pb, i)
    return(sum(zup, zdn, na.rm=T)/sqrt(length(which(!is.na(c(zup, zdn))))))
  }, dset=dset, myeth=myeth, pb=pb)
  tmp[is.na(tmp)] <- 0
  names(tmp) <- rownames(dset)
  return(tmp)
}

myethRea <- function(dset, myeth) {
  myeth[is.na(myeth)] <- 0
  d1 <- t(apply(dset, 1, rank)*2/(ncol(dset)+1)-1)
  es <- rowSums(d1*myeth)/rowSums(abs(myeth))
  return(es*sqrt(rowSums(abs(myeth))))
}



