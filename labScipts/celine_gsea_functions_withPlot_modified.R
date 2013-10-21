gsea2_sep <- function(reflist, gspos, gsneg, np=1000, w=1, sizelim=1, max=F, plot=T,
		pheno1="pos", pheno2="neg", col.f.1="red", col.f.2="blue",
		xlab = "Gene List Index", ylab = "Running Enrichment Score", main.string="GSEA enrichment plots"
) {
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
		gspos <- intersect(names(reflist), as.character(gspos))
		
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
		
		# compute p-value
		if (np>0){
			bg.es.pos = rep(0, np);
			bg.es.pos <- sapply(1:np, function(i){
						bg.isgs <- isgs.pos[sample(1:length(isgs.pos))]  
						bg.hit <- cumsum((abs(reflist*bg.isgs))^w)
						bg.hit <- bg.hit/tail(bg.hit, 1)
						bg.miss <- cumsum(1-bg.isgs)
						bg.miss <- bg.miss/tail(bg.miss, 1)
						bg.all <- bg.hit - bg.miss
						bges <- max(bg.all) + min(bg.all) 
						if(max){
							if(max(bg.all) > -min(bg.all)) bges <- max(bg.all) 
							else bges <- min(bg.all)
						}
						return(bges)    
					})
			if (es.pos<0){
				pv.pos <- length(which(bg.es.pos<=es.pos))/length(which(bg.es.pos<0));
				nes.pos <- es.pos/abs(mean(bg.es.pos[bg.es.pos<0]));
			} else {
				pv.pos <- length(which(bg.es.pos>=es.pos))/length(which(bg.es.pos>0));
				nes.pos <- es.pos/abs(mean(bg.es.pos[bg.es.pos>0]));
			}
		}
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
		
		if (np>0){
			bg.es.neg <- sapply(1:np, function(i){
						bg.isgs <- isgs.neg[sample(1:length(isgs.neg))]
						bg.hit <- cumsum((abs(reflist*bg.isgs))^w)
						bg.hit <- bg.hit/tail(bg.hit, 1)
						bg.miss <- cumsum(1-bg.isgs)
						bg.miss <- bg.miss/tail(bg.miss, 1)
						bg.all <- bg.hit - bg.miss
						bges <- max(bg.all) + min(bg.all) 
						if(max){
							if(max(bg.all) > -min(bg.all)) bges <- max(bg.all) 
							else bges <- min(bg.all)
						}
						return(bges)
					})
			if (es.neg<0){
				pv.neg <- length(which(bg.es.neg<=es.neg))/length(which(bg.es.neg<0))
				nes.neg <- es.neg/abs(mean(bg.es.neg[bg.es.neg<0]))
			} else {
				pv.neg <- length(which(bg.es.neg>=es.neg))/length(which(bg.es.neg>0))
				nes.neg <- es.neg/abs(mean(bg.es.neg[bg.es.neg>0]))
			}
		}
		
	}
	
	if(plot)
	{
		min.RES <- min(min(es_all.pos), min(es_all.neg))
		max.RES <- max(max(es_all.pos),max(es_all.neg))
		if (max.RES < 0.3) max.RES <- 0.3
		if (min.RES > -0.3) min.RES <- -0.3
		#min.RES<- -1.0
		min.RES<- -1.0
		#max.RES<- 0.5
		max.RES<- 1.0
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
		plot(ind, es_all.pos, main = main.string,  xlab = xlab,
				ylab = "Running Enrichment Score", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col.f.1)
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
		
		leg.txt <- pheno1
		text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)
		
		leg.txt <- pheno2
		text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)
		
		adjx.pos <- ifelse(es.pos > 0, 0, 1)
		adjx.neg <- ifelse(es.neg > 0, 0, 1)
		
#    leg.txt <- paste("Peak at ", l.ledge.ref.plot, sep="", collapse="")
		leg.txt <- paste("ES= ", round(es.pos, 2), ", NES= ", round(nes.pos, 2), ", p-value= ", signif(pv.pos, 2), sep="", collapse="")
		text(x=l.ledge.ref.plot.pos, y=max.plot - 0.9*delta, adj = c(adjx.pos, 0), labels=leg.txt, cex = 0.8, col=col.f.1)
		leg.txt <- paste("ES= ", round(es.neg, 2), ", NES= ", round(nes.neg, 2), ", p-value= ", signif(pv.neg, 2), sep="", collapse="")
		text(x=l.ledge.ref.plot.neg, y= min.plot + 1.8*delta, adj = c(adjx.neg, 0), labels=leg.txt, cex = 0.8, col=col.f.2)
		
		#leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
		#text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
		dev.off()
		
		
	}
	
	return(list(es.pos, 
					nes.pos, 
					pv.pos, 
					l.ledge.gs.pos, 
					l.ledge.ref.pos, 
					oddsR.pos,
					es.neg, 
					nes.neg, 
					pv.neg, 
					l.ledge.gs.neg, 
					l.ledge.ref.neg, 
					oddsR.neg,
					ledge.gs.pos, 
					ledge.gs.neg))
}



gsea2_le <- function(reflist, gspos, gsneg, w = 1, sizelim = 1) {
	#GSEA2 Gene set enrichment analysis for 2 complementary gene sets
	#  reflist : named vector of reference scores
	#  gspos   : positive gene set
	#  gsneg   : negative gene set
	#  w     : weight
	#  es    : enrichment score
	#  nes   : normalized enrichment score
	#  pv    : p-value from the permutation test
	#  ledge : leading edge
	#
	# Author: Wei Keat Lim (wl2131@columbia.edu)
	# Modified from Matlab to R by Celine Lefebvre (lefebvre@c2b2.columbia.edu)
	#
	
	#get genes in gs that are in the ref list
	gspos <- intersect(names(reflist), gspos)
	gsneg <- intersect(names(reflist), gsneg)
	
	# combine ranked list and score
	pn <- rep(1, length(reflist))
	ix <- order(reflist, decreasing=T)
	reflist <- reflist[ix]
	pn <- pn[ix]
	
	es.pos <- 0
	l.ledge.gs.pos <- 0
	l.ledge.ref.pos <- 0
	ledge.gs.pos <- NULL
	ledge.ref.pos <- NULL
	
	if(!is.null(gspos) && length(gspos) >= sizelim){
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
#    if(max(es_all.pos) > -min(es_all.pos)) es.pos <- max(es_all.pos) 
#    else es.pos <- min(es_all.pos)
		
		
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
	}      
	
	es.neg <- 0
	l.ledge.gs.neg <- 0
	l.ledge.ref.neg <- 0
	ledge.gs.neg <- NULL
	ledge.ref.neg <- NULL
	
	if(!is.null(gsneg) && length(gsneg) >= sizelim){
		isgs.neg <- rep(0, length(reflist))
		isgs.neg[which(names(reflist) %in% gsneg)] <- 1
		
		score_hit.neg <- cumsum((abs(reflist*isgs.neg))^w)
		score_hit.neg <- score_hit.neg/tail(score_hit.neg, 1)
		score_miss.neg <- cumsum(1-isgs.neg);
		score_miss.neg <- score_miss.neg/tail(score_miss.neg, 1)
		es_all.neg <- score_hit.neg - score_miss.neg
		es.neg <- max(es_all.neg) + min(es_all.neg) 
#    if(max(es_all.neg) > -min(es_all.neg)) es.neg <- max(es_all.neg) 
#    else es.neg <- min(es_all.neg)
		
		
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
			ledge.gs.neg = names(reflist[il1[which(il1 %in% il2)]])
			l.ledge.gs.neg <- length(ledge.gs.neg)
		}
	}
	
	return(list(l.ledge.gs.pos,
					l.ledge.ref.pos,
					l.ledge.gs.neg,
					l.ledge.ref.neg,
					ledge.gs.pos,
					ledge.gs.neg,
					es.pos,
					es.neg))
}

gsea2_es <- function(reflist, gspos, gsneg, w = 1, sizelim = 1) {
	#GSEA2 Gene set enrichment analysis for 2 complementary gene sets
	#  reflist : named vector of reference scores
	#  gspos   : positive gene set
	#  gsneg   : negative gene set
	#  w     : weight
	#  es    : enrichment score
	#
	# Author: Wei Keat Lim (wl2131@columbia.edu)
	# Modified from Matlab to R by Celine Lefebvre (lefebvre@c2b2.columbia.edu)
	#
	
	# combine ranked list and score
	ix <- order(reflist, decreasing=T)
	reflist <- reflist[ix]
	
	es.pos <- 0
	
	if(!is.null(gspos) && length(gspos) >= sizelim)
	{
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
		#if(max(es_all.pos) > -min(es_all.pos)) es.pos <- max(es_all.pos) 
		#else es.pos <- min(es_all.pos)
	}      
	
	es.neg <- 0
	
	if(!is.null(gsneg) && length(gsneg) >= sizelim)
	{
		isgs.neg <- rep(0, length(reflist))
		isgs.neg[which(names(reflist) %in% gsneg)] <- 1
		score_hit.neg <- cumsum((abs(reflist*isgs.neg))^w)
		score_hit.neg <- score_hit.neg/tail(score_hit.neg, 1)
		score_miss.neg <- cumsum(1-isgs.neg);
		score_miss.neg <- score_miss.neg/tail(score_miss.neg, 1)
		es_all.neg <- score_hit.neg - score_miss.neg
		es.neg <- max(es_all.neg) + min(es_all.neg) 
		#if(max(es_all.neg) > -min(es_all.neg)) es.neg <- max(es_all.neg) 
		#else es.neg <- min(es_all.neg)  
	}
	
	return(list(es.pos,
					es.neg))
}

gsea <- function(reflist, gs, np=1000, w=1, sizelim=1, max=F, plot=T, pheno1="", pheno2="", col.f="red",
		xlab = "Gene List Index", ylab = "Running Enrichment Score",main.string="GSEA enrichment plots") {
	#GSEA Gene set enrichment analysis 
	#  reflist : named vector of reference scores
	#  gs   : gene set
	#  np    : number of permutation
	#  w     : weight
	#  es    : enrichment score
	#  nes   : normalized enrichment score
	#  pv    : p-value from the permutation test
	#  ledge : leading edge
	
	
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
			bg.es = rep(0, np);
			bg.es <- sapply(1:np, function(i){
						bg.isgs <- isgs[sample(1:length(isgs))]  
						bg.hit <- cumsum((abs(reflist*bg.isgs))^w)
						bg.hit <- bg.hit/tail(bg.hit, 1)
						bg.miss <- cumsum(1-bg.isgs)
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
		min.RES <- min(es_all)
		max.RES <- max(es_all)
		if (max.RES < 0.3) max.RES <- 0.3
		if (min.RES > -0.3) min.RES <- -0.3
		delta <- (max.RES - min.RES)*0.50
		max.RES<- 0.8
		min.RES<- -0.8
		#max.RES<- 0.5
		#min.RES<- -0.5
		
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
		sub.string <- paste("Number of genes: ", N, " (with measured VIPER Activity), ", length(gs), " (in gene set)", sep = "", collapse="")           
		#main.string <- "Gene Set"
		plot(ind, es_all, main = main.string, sub = sub.string, xlab = xlab, ylab = ylab, xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col.f)
		#for (j in seq(1, N, 20)) {
		#   lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = "grey") # shading of correlation plot
		#}
		lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
		lines(c(l.ledge.ref.plot, l.ledge.ref.plot), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col.f) # max enrichment vertical line
		for (j in 1:N) {
			if (isgs[j]  == 1) {
				lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), 
						lwd = 1, lty = 1, cex = 1, col = "black")  # enrichment tags
			}
		}
		#lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
		lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
		temp <- order(abs(reflist), decreasing=T)
		arg.correl <- temp[N]
		#lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line
		
#		leg.txt <- paste("\"", pheno1, "\" ", sep="", collapse="")
#		text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)
		
#		leg.txt <- paste("\"", pheno2, "\" ", sep="", collapse="")
#		text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)
		
		adjx <- ifelse(es > 0, 0, 1)
		
		#    leg.txt <- paste("Peak at ", l.ledge.ref.plot, sep="", collapse="")
		leg.txt <- paste("NES= ", round(nes, 2), ", p-value= ", signif(pv, 2), sep="", collapse="")
		text(x=l.ledge.ref.plot, y=min.plot + 1.8*delta, adj = c(adjx, 0), 
				labels=leg.txt, cex = 0.8, col="black")
		
		#leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
		#text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
		
		
	}
	return(list(ES=es,
					NES=nes,
					p.value=pv,
					LE.size=l.ledge.gs,
					LE.REF.size=l.ledge.ref,
					LE.genes=ledge.gs,
					oddsRatio=oddsR))
}

gsea.weighted <- function(reflist, gs, np=1000, w=1, sizelim=1, gw) {
	#GSEA Gene set enrichment analysis 
	#  reflist : named vector of reference scores
	#  gs   : gene set
	#  np    : number of permutation
	#  w     : weight
	#  es    : enrichment score
	#  nes   : normalized enrichment score
	#  pv    : p-value from the permutation test
	#  ledge : leading edge
	
	
	#get genes in gs that are in the ref list
	gs <- intersect(names(reflist), gs)
	
	# combine ranked list and score
	pn <- rep(1, length(reflist))
	ix <- order(reflist, decreasing=T)
	reflist <- reflist[ix]
	pn <- pn[ix]
	gw <- 1/gw[ix]
	
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
		isgs.weighted <- isgs
		isgs.weighted[which(names(reflist) %in% gs)] <- gw[which(names(reflist) %in% gs)]
		
		# compute ES
		score_hit <- cumsum((abs(reflist*isgs.weighted))^w)
		score_hit <- score_hit/tail(score_hit, 1)
		score_miss <- cumsum(1-isgs.weighted)
		score_miss <- score_miss/tail(score_miss, 1)
		es_all <- score_hit - score_miss
		es <- max(es_all) + min(es_all) 
		#if(max(es_all) > -min(es_all)) es <- max(es_all) 
		#else es <- min(es_all)
		
		
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
			bg.es = rep(0, np);
			bg.es <- sapply(1:np, function(i){
						perm <- sample(1:length(isgs))
						bg.isgs <- isgs[perm]
						bg.isgs.weighted <- bg.isgs*gw[perm]
						bg.hit <- cumsum((abs(reflist*bg.isgs.weighted))^w)
						bg.hit <- bg.hit/tail(bg.hit, 1)
						bg.miss <- cumsum(1-bg.isgs.weighted)
						bg.miss <- bg.miss/tail(bg.miss, 1)
						bg.all <- bg.hit - bg.miss
						bges <- max(bg.all) + min(bg.all) 
						#if(max(bg.all) > -min(bg.all)) bges <- max(bg.all) 
						#else bges <- min(bg.all)
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
	return(list(es, nes, pv, l.ledge.gs, l.ledge.ref, ledge.gs, oddsR))
}

gsea_es <- function(reflist, gs, w = 1, sizelim = 1) {
	#GSEA Gene set enrichment analysis for 2 complementary gene sets
	#  reflist : named vector of reference scores
	#  gs   : gene set
	#  w     : weight
	#  es    : enrichment score
	#
	# Author: Wei Keat Lim (wl2131@columbia.edu)
	# Modified from Matlab to R by Celine Lefebvre (lefebvre@c2b2.columbia.edu)
	#
	
	# combine ranked list and score
	ix <- order(reflist, decreasing=T)
	reflist <- reflist[ix]
	
	es <- 0
	
	if(!is.null(gs) && length(gs) >= sizelim)
	{
		# check overlap of gene sets and ranked list
		isgs <- rep(0, length(reflist))
		isgs[which(names(reflist) %in% gs)] <- 1
		
		# compute ES
		score_hit <- cumsum((abs(reflist*isgs))^w)
		score_hit <- score_hit/tail(score_hit, 1)
		score_miss <- cumsum(1-isgs)
		score_miss <- score_miss/tail(score_miss, 1)
		es_all <- score_hit - score_miss
		#print(c(max(es_all), min(es_all)))
		es <- max(es_all) + min(es_all) 
		#if(max(es_all) > -min(es_all)) es <- max(es_all) 
		#else es <- min(es_all)
	}  
	return(es)
}

gsea_es_weighted <- function(reflist, gs, w = 1, sizelim = 1, gw) {
	#get genes in gs that are in the ref list
	gs <- intersect(names(reflist), gs)
	
	# combine ranked list and score
	ix <- order(reflist, decreasing=T)
	reflist <- reflist[ix]
	gw <- 1/gw[ix]
	
	es <- 0
	
	if(!is.null(gs) && length(gs) >= sizelim)
	{
		# check overlap of gene sets and ranked list
		isgs <- rep(0, length(reflist))
		isgs[which(names(reflist) %in% gs)] <- 1
		isgs.weighted <- isgs
		isgs.weighted[which(names(reflist) %in% gs)] <- gw[which(names(reflist) %in% gs)]
		
		# compute ES
		score_hit <- cumsum((abs(reflist*isgs))^w)
		score_hit <- score_hit/tail(score_hit, 1)
		score_miss <- cumsum(1-isgs)
		score_miss <- score_miss/tail(score_miss, 1)
		es_all <- score_hit - score_miss
		#print(c(max(es_all), min(es_all)))
		es <- max(es_all) + min(es_all) 
		#if(max(es_all) > -min(es_all)) es <- max(es_all) 
		#else es <- min(es_all)
	}  
	return(es)
}

gsea_le_weighted <- function(reflist, gs, w = 1, sizelim = 1, gw) {
	#get genes in gs that are in the ref list
	gs <- intersect(names(reflist), gs)
	
	# combine ranked list and score
	pn <- rep(1, length(reflist))
	ix <- order(reflist, decreasing=T)
	reflist <- reflist[ix]
	pn <- pn[ix]
	gw <- 1/gw
	
	es <- 0
	l.ledge.gs <- 0
	l.ledge.ref <- 0
	ledge.gs <- NULL
	ledge.ref <- NULL
	
	if(!is.null(gs) && length(gs) >= sizelim){
		# check overlap of gene sets and ranked list 
		isgs <- rep(0, length(reflist))
		isgs[which(names(reflist) %in% gs)] <- 1
		isgs.weighted <- isgs
		isgs.weighted[which(names(reflist) %in% gs)] <- gw[which(names(reflist) %in% gs)]
		
		# compute ES
		score_hit <- cumsum((abs(reflist*isgs.weighted))^w)
		score_hit <- score_hit/tail(score_hit, 1)
		score_miss <- cumsum(1-isgs.weighted)
		score_miss <- score_miss/tail(score_miss, 1)
		es_all <- score_hit - score_miss
		es <- max(es_all) + min(es_all) 
		#if(max(es_all) > -min(es_all)) es <- max(es_all) 
		#else es <- min(es_all)
		
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
	}      
	
	return(list(es,
					l.ledge.gs,
					l.ledge.ref,
					ledge.gs))
}

gsea_le <- function(reflist, gs, w = 1, sizelim = 1) {
	#GSEA Gene set enrichment analysis 
	#  reflist : named vector of reference scores
	#  gspos   : positive gene set
	#  gsneg   : negative gene set
	#  w     : weight
	#  es    : enrichment score
	#  nes   : normalized enrichment score
	#  pv    : p-value from the permutation test
	#  ledge : leading edge
	#
	# Author: Wei Keat Lim (wl2131@columbia.edu)
	# Modified from Matlab to R by Celine Lefebvre (lefebvre@c2b2.columbia.edu)
	#
	
	#get genes in gs that are in the ref list
	gs <- intersect(names(reflist), gs)
	
	# combine ranked list and score
	pn <- rep(1, length(reflist))
	ix <- order(reflist, decreasing=T)
	reflist <- reflist[ix]
	pn <- pn[ix]
	
	es <- 0
	l.ledge.gs <- 0
	l.ledge.ref <- 0
	ledge.gs <- NULL
	ledge.ref <- NULL
	
	if(!is.null(gs) && length(gs) >= sizelim){
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
		#if(max(es_all) > -min(es_all)) es <- max(es_all) 
		#else es <- min(es_all)
		
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
	}      
	
	return(list(es,
					l.ledge.gs,
					l.ledge.ref,
					ledge.gs))
}


getShuffledTtest <- function(mexp, pheno1.or, abs=F, pval=F)
{
	pheno1.length <- length(pheno1.or)
	pheno <- colnames(mexp)
	pheno1 <- sample(pheno, pheno1.length)
	#while(length(intersect(pheno1, pheno1.or))/pheno1.length >= 0.8 || length(intersect(pheno1, pheno1.or))/pheno1.length <= 0.2)
	#  pheno1 <- sample(pheno, pheno1.length)
	
	#str(pheno1)
	pheno2 <- pheno[!(pheno %in% pheno1)]
	#str(pheno2)
	
	if(pval) ttest<-sapply(rownames(mexp), function(p) {
					t<-t.test(mexp[p, pheno1], mexp[p, pheno2]);
					return(sign(t$stat)*(-log(t$p.value, 10)))})
	else ttest <- sapply(rownames(mexp), function(i)
					t.test(mexp[i, pheno1], mexp[i, pheno2])$stat)
	names(ttest) <- rownames(mexp)
	if(abs) ttest <- abs(ttest)
	return(ttest)
}

gsea_pvalueFromES <- function(reflist, gs, bg.es)
{
	es <- 0
	nes <- 0
	pv <- 1
	
	gseale <- gsea_le(reflist, gs)
	es <- gseale[[1]]
	#print(es)
	if(es != 0 )
	{
		# compute p-value
		np = length(bg.es)
		if (np>0){
			if (es<0){
				if(length(bg.es[bg.es<0])>0)
				{
					pv <- length(which((bg.es<=es)))/length(which(bg.es<0));
					nes <- es/mean(abs(bg.es[which(bg.es<0)]))
				}
				else
				{
					pv <- 1
					nes <- 0
					
				}
			}
			else {
				if(length(bg.es[bg.es>0])>0)
				{
					pv <- length(which(bg.es>=es))/length(which(bg.es>0));
					nes <- es/mean(bg.es[which(bg.es>0)])
				}
				else
				{
					pv <- 1
					nes <- 0
				}
			}
		}
	}
	return(list(es,
					nes,
					pv,
					gseale[[2]],
					gseale[[3]]))
}


