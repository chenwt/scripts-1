getNull <- function(mexp=NULL, pheno1.col=NULL, pheno2.col=NULL, reflist=NULL, np=1000, abs=F, pval=F){
  shuffleSamples <- F
  if( (!is.null(mexp)) && ncol(mexp) > 12) {
    shuffleSamples <- T
    pheno <- c(pheno1.col, pheno2.col)
  }
  
  if(shuffleSamples) {
    rnames <- rownames(mexp)
    bg.stat <- sapply(1:np, function(i){
    pheno1.col.random <- sample(pheno, length(pheno1.col))
    pheno2.col.random <- pheno[!pheno %in% pheno1.col.random]
    bg <- myttest(mexp[, pheno2.col.random], mexp[, pheno1.col.random])
    if(pval) bg <- -log(bg$p.value, 10) * sign(bg$stat)
    else bg <- bg$stat
    if(abs) bg <- abs(bg)
    return(bg)
  })}
  else  {
    if(is.null(reflist)) {
      reflist <- myttest(mexp[, pheno2.col], mexp[, pheno1.col])
      if(pval) reflist <- -log(reflist$p.value, 10) * sign(reflist$stat)
      else reflist <- reflist$stat
    }
    if(abs) reflist <- abs(reflist)
    rnames <- names(reflist)
    bg.stat <- sapply(1:np, function(i){
      bg <- reflist[sample(1:length(reflist))]
      return(bg)
    })
  }
  return(bg.stat)
}

gsea2_sep <- function(mexp=NULL, pheno1.col=NULL, pheno2.col=NULL, reflist=NULL, gspos, gsneg, nullDist, w=1, sizelim=1, max=F, plot=T, pheno1="neg", pheno2="pos", col.f.1="red", col.f.2="blue", plotFile="RES.pdf", main.string="GSEA enrichment plots", abs=F, pval=F) {
  #GSEA2 Gene set enrichment analysis for 2 complementary gene sets
  #  reflist : named vector of reference scores
  #  gspos   : positive gene set
  #  gsneg   : negative gene set
  #  np    : number of permutation
  #  w     : weight
  #  es    : enrichment score
  #  nes   : normalized enrichment score
  #  pv    : p-value from the permutation test
  #  ledge : leading edge
  #
  # Author: Wei Keat Lim (wl2131@columbia.edu)
  # Modified from Matlab to R by Celine Lefebvre (lefebvre@c2b2.columbia.edu)
  #

  np <- ncol(nullDist)
  
  if(is.null(mexp) && is.null(reflist)) stop("please provide reference - expression matrix or vector")

  if(is.null(reflist)) {
    reflist <- myttest(mexp[, pheno2.col], mexp[, pheno1.col])
    if(pval) reflist <- -log(reflist$p.value, 10) * sign(reflist$stat)
    else reflist <- reflist$stat
  }
  if(abs) reflist <- abs(reflist)
  
  #get genes in gs that are in the ref list
  gspos <- intersect(names(reflist), gspos)
  gsneg <- intersect(names(reflist), gsneg)

  # combine ranked list and score
  pn <- rep(1, length(reflist))
  ix <- order(reflist, decreasing=T)
  reflist <- reflist[ix]
  pn <- pn[ix]
  
  es.pos <- 0
  nes.pos <- 0 
  pv.pos <- 1
  oddsR.pos <- 0
  l.ledge.gs.pos <- 0
  l.ledge.ref.pos <- 0
  ledge.gs.pos <- NULL
  ledge.ref.pos <- NULL

  if(!is.null(gspos) && length(gspos)>=sizelim){
    # check overlap of gene sets and ranked list 
    isgs.pos <- rep(0, length(reflist))
    isgs.pos[which(names(reflist) %in% gspos)] <- 1
    
    # compute ES
    score_hit.pos <- cumsum((abs(reflist*isgs.pos))^w)
    score_hit.pos <- score_hit.pos/tail(score_hit.pos, 1)
    score_miss.pos <- cumsum(1-isgs.pos)
    score_miss.pos <- score_miss.pos/tail(score_miss.pos, 1)
    es_all.pos <- score_hit.pos - score_miss.pos
    es.pos <- max(es_all.pos) + min(es_all.pos) 
    if(max){
      if(max(es_all.pos) > -min(es_all.pos)) es.pos <- max(es_all.pos) 
      else es.pos <- min(es_all.pos)
    }
    
    # identify leading edge
    isen.pos <- rep(0, length(es_all.pos))
    if (es.pos<0){
        ixpk <- which(es_all.pos==min(es_all.pos))
        isen.pos[ixpk:length(isen.pos)] <- 1
        il1 <- which(isen.pos == 1)
        ledge.ref.pos <- names(reflist[il1])
        l.ledge.ref.pos <- length(ledge.ref.pos)
        il2 <- which(isgs.pos == 1)
        ledge.gs.pos <- names(reflist[il1[which(il1 %in% il2)]])
        ledge.gs.pos <- ledge.gs.pos[length(ledge.gs.pos):1]
        l.ledge.gs.pos <- length(ledge.gs.pos)
    }
    else{
        ixpk <- which(es_all.pos==max(es_all.pos))
        isen.pos[1:ixpk] <- 1
        il1 <- which(isen.pos == 1)
        ledge.ref.pos <- names(reflist[il1])
        l.ledge.ref.pos <- length(ledge.ref.pos)
        il2 <- which(isgs.pos == 1)
        ledge.gs.pos = names(reflist[il1[which(il1 %in% il2)]])
        l.ledge.gs.pos <- length(ledge.gs.pos)
    }
    
    oddsR.pos <- (l.ledge.gs.pos/l.ledge.ref.pos)/((length(gspos)-l.ledge.gs.pos)/(length(reflist)-l.ledge.ref.pos))
  }
  
  es.neg <- 0
  nes.neg <- 0 
  pv.neg <- 1
  oddsR.neg <- 1
  l.ledge.gs.neg <- 0
  l.ledge.ref.neg <- 0
  ledge.gs.neg <- NULL
  ledge.ref.neg <- NULL
  
  if(!is.null(gsneg) && length(gsneg)>=sizelim){
    gsneg <- intersect(names(reflist), as.character(gsneg))
    isgs.neg <- rep(0, length(reflist))
    isgs.neg[which(names(reflist) %in% gsneg)] <- 1
    
    score_hit.neg <- cumsum((abs(reflist*isgs.neg))^w)
    score_hit.neg <- score_hit.neg/tail(score_hit.neg, 1)
    score_miss.neg <- cumsum(1-isgs.neg);
    score_miss.neg <- score_miss.neg/tail(score_miss.neg, 1)
    es_all.neg <- score_hit.neg - score_miss.neg
    es.neg <- max(es_all.neg) + min(es_all.neg) 
    if(max){
      if(max(es_all.neg) > -min(es_all.neg)) es.neg <- max(es_all.neg) 
      else es.neg <- min(es_all.neg)
    }
    
    isen.neg <- rep(0, length(es_all.neg))
    if (es.neg<0){
      ixpk <- which(es_all.neg==min(es_all.neg))
      isen.neg[ixpk:length(isen.neg)] <- 1
      il1 <- which(isen.neg == 1)
      ledge.ref.neg <- names(reflist[il1])
      l.ledge.ref.neg <- length(ledge.ref.neg)
      il2 <- which(isgs.neg == 1)
      ledge.gs.neg <- names(reflist[il1[which(il1 %in% il2)]])
      ledge.gs.neg <- ledge.gs.neg[length(ledge.gs.neg):1]
      l.ledge.gs.neg <- length(ledge.gs.neg)
    }
    else{
      ixpk <- which(es_all.neg==max(es_all.neg))
      isen.neg[1:ixpk] = 1
      il1 <- which(isen.neg == 1)
      ledge.ref.neg <- names(reflist[il1])
      l.ledge.ref.neg <- length(ledge.ref.neg)
      il2 <- which(isgs.neg == 1)
      ledge.gs.neg <- names(reflist[il1[which(il1 %in% il2)]])
      l.ledge.gs.neg <- length(ledge.gs.neg)
    }
    oddsR.neg <- (l.ledge.gs.neg/l.ledge.ref.neg)/((length(gsneg)-l.ledge.gs.neg)/(length(reflist)-l.ledge.ref.neg))
  }
  
  # compute p-values
  if (np>0){
    bg.es <- sapply(1:ncol(nullDist), function(i){
      reflist.random <- nullDist[,i]
      names(reflist.random) <- rownames(nullDist)
      reflist.random <- sort(reflist.random, decreasing=T)

      isgs.pos.random <- rep(0, length(reflist.random))
      isgs.pos.random[which(names(reflist.random) %in% gspos)] <- 1
      
      bg.hit.p <- cumsum((abs(reflist.random*isgs.pos.random))^w)
      bg.hit.p <- bg.hit.p/tail(bg.hit.p, 1)
      bg.miss.p <- cumsum(1-isgs.pos)
      bg.miss.p <- bg.miss.p/tail(bg.miss.p, 1)
      bg.all.p <- bg.hit.p - bg.miss.p
      bges.p <- max(bg.all.p) + min(bg.all.p) 
      if(max){
        if(max(bg.all.p) > -min(bg.all.p)) bges.p <- max(bg.all.p) 
        else bges.p <- min(bg.all.p)
      }
      
      isgs.neg.random <- rep(0, length(reflist.random))
      isgs.neg.random[which(names(reflist.random) %in% gsneg)] <- 1
      
      bg.hit.n <- cumsum((abs(reflist.random*isgs.neg.random))^w)
      bg.hit.n <- bg.hit.n/tail(bg.hit.n, 1)
      bg.miss.n <- cumsum(1-isgs.neg)
      bg.miss.n <- bg.miss.n/tail(bg.miss.n, 1)
      bg.all.n <- bg.hit.n - bg.miss.n
      bges.n <- max(bg.all.n) + min(bg.all.n) 
      if(max){
        if(max(bg.all.n) > -min(bg.all.n)) bges.n <- max(bg.all.n) 
        else bges.n <- min(bg.all.n)
      }
      return(c(bges.p, bges.n))    
    })
    
    if (es.pos<0){
      pv.pos <- length(which(bg.es[1,]<=es.pos))/length(which(bg.es[1,]<0));
      nes.pos <- es.pos/abs(mean(bg.es[1,][bg.es[1,]<0]));
    } else {
      pv.pos <- length(which(bg.es[1,]>=es.pos))/length(which(bg.es[1,]>0));
      nes.pos <- es.pos/abs(mean(bg.es[1,][bg.es[1,]>0]));
    }
    
    if (es.neg<0){
      pv.neg <- length(which(bg.es[2,]<=es.neg))/length(which(bg.es[2,]<0))
      nes.neg <- es.neg/abs(mean(bg.es[2,][bg.es[2,]<0]))
    } else {
      pv.neg <- length(which(bg.es[2,]>=es.neg))/length(which(bg.es[2,]>0))
      nes.neg <- es.neg/abs(mean(bg.es[2,][bg.es[2,]>0]))
    }
  }
  
  if(pv.pos==0) pv.pos <- 1/np
  if(pv.neg==0) pv.neg <- 1/np
  

  if(plot)
  {
    pdf(plotFile)
    min.RES <- min(min(es_all.pos), min(es_all.neg))
    max.RES <- max(max(es_all.pos),max(es_all.neg))
    if (max.RES < 0.3) max.RES <- 0.3
    if (min.RES > -0.3) min.RES <- -0.3
    delta <- (max.RES - min.RES)*0.50
    min.plot <- min.RES - 2*delta
    max.plot <- max.RES + delta
    max.corr <- max(reflist)
    min.corr <- min(reflist)
    Obs.correl.vector.norm <- (reflist - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
    zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
    
    if(es.pos<0) l.ledge.ref.plot.pos <- length(reflist)-l.ledge.ref.pos
    else l.ledge.ref.plot.pos <- l.ledge.ref.pos
    if(es.neg<0) l.ledge.ref.plot.neg <- length(reflist)-l.ledge.ref.neg
    else l.ledge.ref.plot.neg <- l.ledge.ref.neg
    
    N<-length(reflist)
    ind <- 1:N
    #sub.string <- paste("Number of genes: ", N, " (in list), ", length(gs), " (in gene set)", sep = "", collapse="")           
    #main.string <- "Gene Set"
    plot(ind, es_all.pos, main = main.string,  xlab = "Gene List Index", ylab = "Running Enrichment Score", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col.f.1)
    points(ind, es_all.neg, type = "l", lwd = 2, cex = 1, col = col.f.2)
    for (j in seq(1, N, 20)) {
      lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = "lightgrey") # shading of correlation plot
    }
    lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
    lines(c(l.ledge.ref.plot.pos, l.ledge.ref.plot.pos), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col.f.1) # max enrichment vertical line
    lines(c(l.ledge.ref.plot.neg, l.ledge.ref.plot.neg), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col.f.2) # max enrichment vertical line
    for (j in 1:N) {
      if (isgs.pos[j]  == 1) {
        lines(c(j, j), c(max.plot - 0.75*delta, max.plot - 0.15*delta), lwd = 1, lty = 1, cex = 1, col = col.f.1)  # enrichment tags
      }
      if (isgs.neg[j]  == 1) {
        lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = col.f.2)  # enrichment tags
      }
    }
    lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
    lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
    temp <- order(abs(reflist), decreasing=T)
    arg.correl <- temp[N]
    #lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line
    
    leg.txt <- pheno2
    text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)
    
    leg.txt <- pheno1
    text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)
    
    adjx.pos <- ifelse(es.pos > 0, 0, 1)
    adjx.neg <- ifelse(es.neg > 0, 0, 1)
           
#    leg.txt <- paste("Peak at ", l.ledge.ref.plot, sep="", collapse="")
        if(pv.pos==0) pv.pos <- 1/np
#    leg.txt <- paste("ES= ", round(es.pos, 2), ", NES= ", round(nes.pos, 2), ", p-value= ", signif(pv.pos, 2), sep="", collapse="")
    leg.txt <- paste("NES= ", round(nes.pos, 2), ", p-value= ", signif(pv.pos, 2), sep="", collapse="")
    text(x=l.ledge.ref.plot.pos, y=max.plot - 0.9*delta, adj = c(adjx.pos, 0), labels=leg.txt, cex = 0.8, col=col.f.1)
#    leg.txt <- paste("ES= ", round(es.neg, 2), ", NES= ", round(nes.neg, 2), ", p-value= ", signif(pv.neg, 2), sep="", collapse="")
        if(pv.neg==0) pv.neg <- 1/np
    leg.txt <- paste("NES= ", round(nes.neg, 2), ", p-value= ", signif(pv.neg, 2), sep="", collapse="")
    text(x=l.ledge.ref.plot.neg, y= min.plot + 1.8*delta, adj = c(adjx.neg, 0), labels=leg.txt, cex = 0.8, col=col.f.2)
    
    #leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
    #text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
    dev.off()
    
    
  }

  return(list(ES.pos=es.pos, 
              NES.pos=nes.pos, 
              p.value.pos=pv.pos, 
              LE.size.pos=l.ledge.gs.pos, 
              LE.ref.size.pos=l.ledge.ref.pos, 
              oddsR.pos=oddsR.pos,
              ES.neg=es.neg, 
              NES.neg=nes.neg, 
              p.value.neg=pv.neg, 
              LE.size.neg=l.ledge.gs.neg, 
              LE.ref.size.neg=l.ledge.ref.neg, 
              oddsR.neg=oddsR.neg,
              LE.genes.pos=ledge.gs.pos, 
              LE.genes.neg=ledge.gs.neg))
}

gsea <- function(mexp=NULL, pheno1.col=NULL, pheno2.col=NULL, reflist=NULL,  gs, nullDist, w=1, sizelim=1, max=F, plot=T, pheno1="neg", pheno2="pos", col.f="red", plotFile="RES.pdf", main.string="GSEA enrichment plots", abs=F, pval=F) {
  #GSEA Gene set enrichment analysis 
  #  reflist : named vector of reference scores
  #  gs   : gene set
  #  np    : number of permutation
  #  w     : weight
  #  es    : enrichment score
  #  nes   : normalized enrichment score
  #  pv    : p-value from the permutation test
  #  ledge : leading edge

  np <- ncol(nullDist)
  
  if(is.null(mexp) && is.null(reflist)) stop("please provide reference - expression matrix or vector")

  if(is.null(reflist)) {
    reflist <- myttest(mexp[, pheno2.col], mexp[, pheno1.col])
    if(pval) reflist <- -log(reflist$p.value, 10) * sign(reflist$stat)
    else reflist <- reflist$stat
  }
  if(abs) reflist <- abs(reflist)
  
  if((length(pheno1.col)+length(pheno2.col))>12) shuffleSamples <- T

  #get genes in gs that are in the ref list
  gs <- intersect(names(reflist), gs)
  
  # combine ranked list and score
  pn <- rep(1, length(reflist))
  ix <- order(reflist, decreasing=T)
  reflist <- reflist[ix]
  pn <- pn[ix]
  
  es <- 0
  nes <- 0 
  pv <- 1
  l.ledge.gs <- 0
  l.ledge.ref <- 0
  ledge.gs <- NULL
  ledge.ref <- NULL
  oddsR <- 0

  if(!is.null(gs) && length(gs)>=sizelim){
    # check overlap of gene sets and ranked list 
    isgs <- rep(0, length(reflist))
    isgs[which(names(reflist) %in% gs)] <- 1
      
    # compute ES
    score_hit <- cumsum((abs(reflist*isgs))^w)
    score_hit <- score_hit/tail(score_hit, 1)
    score_miss <- cumsum(1-isgs)
    score_miss <- score_miss/tail(score_miss, 1)
    es_all <- score_hit - score_miss
    es <- max(es_all) + min(es_all) 
    if(max){
      if(max(es_all) > -min(es_all)) es <- max(es_all) 
      else es <- min(es_all)
    }
    
    # identify leading edge
    isen <- rep(0, length(es_all))
    if (es<0){
        ixpk <- which(es_all==min(es_all))
        isen[ixpk:length(isen)] <- 1
        il1 <- which(isen == 1)
        ledge.ref <- names(reflist[il1])
        l.ledge.ref <- length(ledge.ref)
        il2 <- which(isgs == 1)
        ledge.gs <- names(reflist[il1[which(il1 %in% il2)]])
        ledge.gs <- ledge.gs[length(ledge.gs):1]
        l.ledge.gs <- length(ledge.gs)
    }
    else{
        ixpk <- which(es_all==max(es_all))
        isen[1:ixpk] <- 1
        il1 <- which(isen == 1)
        ledge.ref <- names(reflist[il1])
        l.ledge.ref <- length(ledge.ref)
        il2 <- which(isgs == 1)
        ledge.gs = names(reflist[il1[which(il1 %in% il2)]])
        l.ledge.gs <- length(ledge.gs)
    }
      
    # compute p-value
    if (np>0){
      bg.es <- sapply(1:ncol(nullDist), function(i){
        reflist.random <- nullDist[,i]
        names(reflist.random) <- rownames(nullDist)
        reflist.random <- sort(reflist.random, decreasing=T)
        isgs.random <- rep(0, length(reflist.random))
        isgs.random[which(names(reflist.random) %in% gs)] <- 1

        bg.hit <- cumsum((abs(reflist.random*isgs.random))^w)
        bg.hit <- bg.hit/tail(bg.hit, 1)
        bg.miss <- cumsum(1-isgs)
        bg.miss <- bg.miss/tail(bg.miss, 1)
        bg.all <- bg.hit - bg.miss
        bges <- max(bg.all) + min(bg.all) 
        if(max){
          if(max(bg.all) > -min(bg.all)) bges <- max(bg.all) 
          else bges <- min(bg.all)
        }
        return(bges)    
      })
      
      if (es<0){
        pv <- sum(bg.es<=es)/length(which(bg.es<0));
        nes <- es/abs(mean(bg.es[bg.es<0]));
      } else {
        pv <- sum(bg.es>=es)/length(which(bg.es>0));
        nes <- es/abs(mean(bg.es[bg.es>0]));
      }
    }
    oddsR <- (l.ledge.gs/l.ledge.ref)/((length(gs)-l.ledge.gs)/(length(reflist)-l.ledge.ref))
  }  

  if(pv==0) pv <- 1/np
  
  if(plot)
  {
    pdf(plotFile)
    min.RES <- min(es_all)
    max.RES <- max(es_all)
    if (max.RES < 0.3) max.RES <- 0.3
    if (min.RES > -0.3) min.RES <- -0.3
    delta <- (max.RES - min.RES)*0.50
    min.plot <- min.RES - 2*delta
    max.plot <- max.RES
    max.corr <- max(reflist)
    min.corr <- min(reflist)
    Obs.correl.vector.norm <- (reflist - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
    zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
    
    if(es<0) l.ledge.ref.plot <- length(reflist)-l.ledge.ref
    else l.ledge.ref.plot <- l.ledge.ref
    

    if(nes>0) col.f <- "red"
    else col.f <- "blue"
    N<-length(reflist)
    ind <- 1:N
    sub.string <- paste("Number of genes: ", N, " (in list), ", length(gs), " (in gene set)", sep = "", collapse="")           
    #main.string <- "Gene Set"
    plot(ind, es_all, main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col.f)
    for (j in seq(1, N, 20)) {
      lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = "grey") # shading of correlation plot
    }
    lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
    lines(c(l.ledge.ref.plot, l.ledge.ref.plot), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col.f) # max enrichment vertical line
    for (j in 1:N) {
      if (isgs[j]  == 1) {
        lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = col.f)  # enrichment tags
      }
    }
    lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
    lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
    temp <- order(abs(reflist), decreasing=T)
    arg.correl <- temp[N]
    #lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line
    
    leg.txt <- paste("\"", pheno2, "\" ", sep="", collapse="")
    text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)
    
    leg.txt <- paste("\"", pheno1, "\" ", sep="", collapse="")
    text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)
    
    adjx <- ifelse(es > 0, 0, 1)
           
#    leg.txt <- paste("Peak at ", l.ledge.ref.plot, sep="", collapse="")
    leg.txt <- paste("NES= ", round(nes, 2), ", p-value= ", signif(pv, 2), sep="", collapse="")
    text(x=l.ledge.ref.plot, y=min.plot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 0.8, col=col.f)
    
    #leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
    #text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
    dev.off()
    
    
  }
  return(list(ES=es,
              NES=nes,
              p.value=pv,
              LE.size=l.ledge.gs,
              LE.REF.size=l.ledge.ref,
              LE.genes=ledge.gs,
              oddsRatio=oddsR))
}

myttest <- function (x, y,  alternative = "two.sided", welch=T)
{
    lx <- ncol(x)
    ly <- ncol(y)
    x.var <- rowVars(x)
    y.var <- rowVars(y)
    if(welch) {
      t <- as.vector((rowMeans(x) - rowMeans(y))/sqrt((x.var/lx+y.var/ly)))
      df <- as.vector((x.var/lx+y.var/ly)^2/( (x.var/lx)^2/(lx-1) + (y.var/ly)^2/(ly-1)))
      df[which(df >  (lx+ly-2))] <- lx+ly-2
      df[which(df <  min(lx,ly))] <- min(lx, ly)
    }
    else{
      t <- as.vector((rowMeans(x) - rowMeans(y))/
        (sqrt(((lx - 1) * x.var + (ly - 1) * y.var)/(lx + ly - 2))*sqrt(1/lx + 1/ly)))
      df <- lx + ly -2
    }
    names(t) <- rownames(x)
    p <- as.vector(switch(pmatch(alternative,
           c("two.sided", "greater", "less")),
           pt(abs(t), df, lower.tail = F) * 2,
           pt(t, df, lower.tail = F),
           pt(t, df, lower.tail = T)))
    names(p) <- rownames(x)
    list(statistic = t,
         p.value = p)
}

rowVars <- function (x)
{
  rowSums((x - rowMeans(x))^2)/ncol(x)
}

