myDunnIndex <- function(dataMA,k=3,dmethod="euclidean",aggmethod="min"){
    mod <- kmeans(dataMA,k)
    cCenters <- mod$centers # matrix for all centroied
    cSize <- mod$size #integer
    cDist <- matrix(NA,nrow = k,ncol = k)
    cMax <- list()  # matrix for each cluster

    ## extract the matrix for each cluster
    for (i in 1:k){
        cMax[[i]] <- dataMA[names(which(mod$cluster == i)),] 
    }
    ## computing inter-cluster distance
    for ( i in 1:k){
        for( j in i:k){
            if(aggmethod == "centroid") {
                cDist[i,j] <- mydist(cCenters[i,],cCenters[j,],dmethod)
            } else {
                cDist[i,j] <- myCdist(cMax[[i]],cMax[[j]],dmethod,aggmethod)
            }
        }
    }

    #### computing set diameter
    delta <- NA
    for (i in 1:k){
        # mean distance to cluster center is 
        # delta[i] <- sum(mydist(cMax[[i]],cCenters[i,],dmethod)) / cSize[i]
        delta[i] <- myCdist(cMax[[i]],cMax[[i]],dmethod,"max") 
    }
    ### computing dunn index
    di2 <- min(cDist[!is.na(cDist)])/max(delta)                    
    return(di2)
}


## this is to calculate the distance(euclidean) between two clusters
myCdist <- function(x,y,dmethod="euclidean",aggmethod="min"){

    nr <- nrow(x)
    nc <- nrow(y)

    if(is.null(nc)) nc <- 1
    if(is.null(nr)) nr <- 1
    cd <- NA
    if(nc > 1 && nr > 1){
        for ( i in (1:nc)){
            switch(aggmethod,
                min = {switch(dmethod,
                          euclidean = { tmpd <- min(mydist(x,y[i,],"euclidean")); if(i == 1) {cd = tmpd}; cd <- min(cd, tmpd) },
                          spearman = { tmpd <- min(mydist(x,y[i,],"spearman")); if(i == 1) {cd = tmpd}; cd <- min(cd, tmpd)  },
                          stop("please enter distance method: euclidean/spearman")
                          )},
                max = {switch(dmethod,
                          euclidean = {tmpd <- max(mydist(x,y[i,],"euclidean"));if(i==1){cd = tmpd}; cd <- max(cd,tmpd)},
                          spearman = {tmpd <- max(mydist(x,y[i,],"spearman")); if(i==1){cd= tmpd};cd <- max(cd,tmpd)},
                          stop("please enter distance method: euclidean/spearman")
                        )
                    },
                stop("please enter inter cluster methods:min/max")
                )
        }
    } else {
        if( nc/nr > 1) {tmp <- x;x <- y; y <- tmp}
            switch(aggmethod,
            min={switch(dmethod,
                    euclidean = {cd <- min(mydist(x,y,"euclidean"))},
                    spearman = {cd <- min(mydist(x,y,"spearman"))},
                    )},
            max={switch(dmethod,
                    euclidean = {cd <- max(mydist(x,y,"euclidean"))},
                    spearman = {cd <- max(mydist(x,y,"spearman"))},
                    )},
            )
    }
    return(cd)
}

### this is to calculated the global distance of each points to centroids( n to 1)
mydist <- function(x,y,method){
    nr <- nrow(x)
    nn <- length(y)
    if(is.null(nr)) nr <- 1
    # for ( i in (1:nc)){
    if(nr > 1) {
        d <- matrix(NA,nrow=nr,ncol=1)
        switch(method,
            euclidean = {d[,1] <- sqrt(rowSums((myMaOps(x,y,"mi"))^2))},
            spearman = {
                xr <- t(apply(x,1,FUN=rank));
                # yr <- t(apply(y,1,FUN=rank));
                # d[,i] <- 1 - (6 * rowSums(myMaOps(xr,rank(y),"mi")^2) / (nn * ( nn ^2 - 1)))},
                d[,1] <- rowSums(myMaOps(xr,rank(y),"mi")^2)},
                # d[,i] <- 1 - (6 * (rowSums((xr-rank(y))^2)) / (nn * ( nn ^2 - 1)))},
            stop("Please enter euclidean or spearman!")
        )
    }else{
        switch(method,
            euclidean = {d <- sqrt(sum((x-y)^2))},
            # spearman = {d <- 1 - (6 * sum((rank(x)-rank(y))^2) /( nn * ( nn ^2 - 1)))},
            spearman = {d <- sum((rank(x)-rank(y))^2)},
            stop("Please enter euclidean or spearman!")
        )
    }
    return(d)
}

myMaOps <- function(x,y,ops) {
    result <- NULL
    switch(ops,
     p = {result <- t(apply(x,1,function(x){(x+y)}))},
     mi = {result <- t(apply(x,1,function(x){(x-y)}))},
     ml = {result <- t(apply(x,1,function(x){(x*y)}))},
     d= {result <- t(apply(x,1,function(x){(x/y)}))},
     stop("Please enter p/+,mi/-,ml/*,d/")
     )
    return (result)
}


myDiBatch <- function(dataMA,mink,maxk,dmethod,aggmethod){
    diDF <- data.frame(cbind(k=mink:maxk,di=NA),dmethod=rep(paste(dmethod),(maxk-mink)),deltamethod=rep(paste(aggmethod),(maxk-mink+1)))

    for(i in mink:maxk)
     diDF[(i-mink+1),2]<- as.numeric(myDunnIndex(dataMA,i,dmethod,aggmethod))
    return(diDF)
}

diDF <- myDiBatch(exp177MA,2,60,"euclidean","max")
diDF <- rbind(diDF,myDiBatch(exp177MA,2,60,"euclidean","min"))
diDF <- rbind(diDF,myDiBatch(exp177MA,2,60,"euclidean","centroid"))
diDF <- rbind(diDF,myDiBatch(exp177MA,2,60,"spearman","max"))
diDF <- rbind(diDF,myDiBatch(exp177MA,2,60,"spearman","min"))

diNew <- cbind(k = seq(2,60),spearman_max = diDF[which(diDF$dmethod == "spearman" & diDF$deltamethod == "max"),2])
diNew <- cbind(diNew,euclidean_max = diDF[which(diDF$dmethod == "euclidean" & diDF$deltamethod == "max"),2])
diNew <- cbind(diNew,euclidean_min = diDF[which(diDF$dmethod == "euclidean" & diDF$deltamethod == "min"),2])



save(diDF,exp.data,exp177MA,exp13MA,file="DunnIndex.RData")

diDF <- rbind(diDF,myDiBatch(exp177MA,2,60,"spearman","centroid"))

ggplot(diDF,aes(x=k,y=di,shape=dmethod,color=deltamethod)) + geom_point(size=2)+geom_line() + scale_x_discrete(breaks = 2:60,labels = 2:60)





ggplot(diDF,aes(x=k, y=di,color=c(dmethod,deltamethod))) + geom_line() + geom_point(size=3)

/ifs/home/c2b2/ac_lab/jh3283/schome/AML/r_process/DunnIndex.RData
