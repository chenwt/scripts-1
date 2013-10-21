myPamrGene<- function (fit, data, threshold) 
{
    par(pch = 1, col = 1)
    geneid <- data$geneid
    if (is.null(geneid)) {
        geneid <- as.character(1:nrow(data$x))
    }
    if (is.null(fit$newy)) {
        y <- factor(data$y[fit$sample.subset])
    }
    else {
        y <- factor(fit$newy[fit$sample.subset])
    }
    x <- data$x[fit$gene.subset, fit$sample.subset]
    geneid <- geneid[fit$gene.subset]
    nc <- length(unique(y))
    aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
    cen <- pamr.predict(fit, x, threshold = threshold, type = "cen")
    
    d <- (cen - fit$centroid.overall)[aa, ]/fit$sd[aa]
    oo <- order(-apply(abs(d), 1, max))
    aa <- aa[oo]
    ngenes <- length(aa)
    o <- order(y)
    xx <- x[aa, o]
    geneid <- geneid[aa]
    
    nc <- length(unique(y))
    nn <- c(0, cumsum(table(y)))
    nrow <- trunc(sqrt(ngenes)) + 1
    ncol <- trunc(sqrt(ngenes)) + 1
    return(list(xx=xx,genelist=as.character(geneid)))
}


myPamrConfusion <- function (fit, threshold, extra = TRUE) 
{
    ii <- (1:length(fit$threshold))[fit$threshold >= threshold]
    ii <- ii[1]
    predicted <- fit$yhat[, ii]
    if (!is.null(fit$y)) {
        true <- fit$y[fit$sample.subset]
        tt <- table(true, predicted)
    }
    else {
        true <- fit$proby[fit$sample.subset, ]
        ytemp <- apply(true, 1, which.is.max)
        temp <- c(predicted, names(table(ytemp)))
        nams <- names(table(temp))
        Yhat <- model.matrix(~factor(temp) - 1, data = list(y = temp))
        Yhat <- Yhat[1:length(predicted), ]
        tt <- matrix(NA, nrow = length(fit$prior), ncol = length(fit$prior))
        for (i in 1:length(fit$prior)) {
            for (j in 1:length(fit$prior)) {
                tt[i, j] <- sum(true[, i] * Yhat[, j])
            }
        }
        dimnames(tt) <- list(names(table(ytemp)), nams)
    }
    if (extra) {
        tt1 <- tt
        diag(tt1) <- 0
        tt <- cbind(tt, apply(tt1, 1, sum)/apply(tt, 1, sum))
        dimnames(tt)[[2]][ncol(tt)] <- "Class Error rate"
        return(tt)
    }
    if (!extra) {
        return(tt)
    }
}


myPlotCV <- function (fit) 
{
    par(mar = c(5, 5, 5, 5),xpd=TRUE)
    par(mfrow = c(2, 1))
    n <- nrow(fit$yhat)
    y <- fit$y
    if (!is.null(fit$newy)) {
        y <- fit$newy[fit$sample.subset]
    }
    nc <- length(table(y))
    nfolds <- length(fit$folds)
    err <- matrix(NA, ncol = ncol(fit$yhat), nrow = nfolds)
    temp <- matrix(y, ncol = ncol(fit$yhat), nrow = n)
    ni <- rep(NA, nfolds)
    for (i in 1:nfolds) {
        ii <- fit$folds[[i]]
        ni[i] <- length(fit$folds[[i]])
        err[i, ] <- apply(temp[ii, ] != fit$yhat[ii, ], 2, sum)/ni[i]
    }
    se <- sqrt(apply(err, 2, var)/nfolds)
    plot(fit$threshold, fit$error, ylim = c(-0.1, 0.8), xlab = "Value of threshold  ", ylab = "Misclassification Error", type = "n", yaxt = "n",col="blue")
    axis(3, at = fit$threshold, lab = paste(fit$size), srt = 90, adj = 0)
    mtext("Number of genes", 3, 4, cex = 1.2)
    axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
    lines(fit$threshold, fit$error, col = 2)
    o <- fit$err == min(fit$err)
    points(fit$threshold[o], fit$error[o], pch = "x")
    #error.bars(fit$threshold, fit$err - se, fit$err + se)
    err2 <- matrix(NA, nrow = length(unique(y)), ncol = length(fit$threshold))
    for (i in 1:(length(fit$threshold) - 1)) {
        s <- pamr.confusion(fit, fit$threshold[i], extra = FALSE)
        diag(s) <- 0
        err2[, i] <- apply(s, 1, sum)/table(y)
    }
    plot(fit$threshold, err2[1, ], ylim = c(-0.1, 1.1), xlab = "Value of threshold ", ylab = "Misclassification Error", type = "n", yaxt = "n")
    axis(3, at = fit$threshold, lab = paste(fit$size), srt = 90, adj = 0)
    axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
    for (i in 1:nrow(err2)) {
        lines(fit$threshold, err2[i, ], col = i + 1)
    }
    legend("right",inset=c(-0.2,0),dimnames(table(y))[[1]], col = (2:(nc + 1)), lty = 1)
    par(mfrow = c(1, 1))
}

source("C:/Documents and Settings/jh3283/My Documents/Dropbox/Viz/heatmap/heatmap3.R")


get_heatmapMX <- function(geneset,data) {
    geneset <- unlist(geneset)
    row.names(data) <- dataPM$genenames
    data.heatmap <- as.matrix(data[geneset,])
    fab <- clinic.data$fab[match(colnames(data.heatmap),clinic.data$usi)]
    uniqfab <- names(table(fab))
    
    event <- clinic.data$efs.event.type.id[match(colnames(data.heatmap),clinic.data$usi)]
    uniqevent <- names(table(event))
    colors <- c("blue","green","yellow","orange","red","purple","black","pink","gray") 
        fab <- gsub("M0",colors[1],fab)
        fab <- gsub("M1",colors[2],fab)
        fab <- gsub("M2",colors[3],fab)
        fab <- gsub("M4",colors[4],fab)
        fab <- gsub("M5",colors[5],fab)
        fab <- gsub("M6",colors[6],fab)
        fab <- gsub("M7",colors[7],fab)
        fab <- gsub("NOS",colors[9],fab)
        fab <- gsub("Unknown",colors[9],fab)

        event <- gsub("Censored",colors[1],event)
        event <- gsub("Death without remission",colors[2],event)
        event <- gsub("Induction failure",colors[3],event)
        event <- gsub("Relapse",colors[4],event)
        event <- gsub("Death",colors[5],event)
    
    return(list(data=data.heatmap,col.fab=fab,fab=uniqfab,col.event=event,event=uniqevent))
}