## this R script is to using Dunn Index to evalute cluster Robustness after doing k-means
myDunnIndex2 <- function(dataMA,k=3,aggmethod="min"){
    mod <- kmeans(dataMA,k)
    cSize <- mod$size #integer
    cDist <- matrix(NA,nrow = k,ncol = k)
    cMax <- list()  # matrix for each cluster
## extract the matrix for each cluster
    for (i in 1:k){
        cMax[[i]] <- t(apply(dataMA[names(which(mod$cluster == i)),],1,FUN=rank))
    }
 ## computing inter-cluster distance
    for ( i in 1:k){
        for( j in 1:k){
                if(i != j ){ cDist[i,j] <- mySCdist(cMax[[i]],cMax[[j]],aggmethod)}
            }
}
#### computing set diameter
    delta <- NA
    for (i in 1:k){
        delta[i] <- mySCdist(cMax[[i]],cMax[[i]],"max") 
    }
### computing dunn index
    di2 <- min(cDist[!is.na(cDist)])/max(delta)                    
    return(di2)
}


##### this is to calculate spearman distance between two matrix
mySCdist <- function(x,y,aggmethod="min"){
    nr <- nrow(x)
    nc <- nrow(y)
    if(is.null(nc)) nc <- 1
    if(is.null(nr)) nr <- 1
    cd <- NA
    if(nc > 1 && nr > 1){
        for ( i in (1:nc)){
            switch(aggmethod,
                min = {tmpd <- min(mySdist(x,y[i,])); if(i == 1) {cd = tmpd}; cd <- min(cd, tmpd)  },
                max = {tmpd <- max(mySdist(x,y[i,])); if(i == 1){cd= tmpd};cd <- max(cd,tmpd)},
                stop("please enter inter cluster methods:min/max")
                )
        }
    } else {
        if( nc/nr > 1) {tmp <- x;x <- y; y <- tmp}
            switch(aggmethod,
            min={cd <- min(mySdist(x,y))},
            max={cd <- max(mySdist(x,y))},
            )
    }
    return(cd)
}

### calculate the spearman distance to one centroid y
mySdist <- function(x,y){
    nr <- nrow(x)
    nn <- length(y)
    if(is.null(nr)) nr <- 1
    if(nr > 1) {
        d <- matrix(NA,nrow=nr,ncol=1)
            d[,1] <- rowSums(myMaOps(x,y,"mi")^2)
    }else{
        d <- sum(x-y)^2
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


myDiBatch2 <- function(dataMA,mink,maxk,aggmethod){
    diDF <- data.frame(cbind(k=mink:maxk,di=NA))
    for(i in mink:maxk)
     diDF[(i-mink+1),2]<- as.numeric(myDunnIndex2(dataMA,i,aggmethod))
    return(diDF)
}



diDF_S <- myDiBatch(exp177MA,2,60,"min")
ggplot(diDF_s,aes(x=k,y=di)) + geom_line(color="blue") + opt(title = "DI", subtitle = "spearman")

diDF_A <- rbind(cbind(diDF,dist = rep("Euclidean",nrow(diDF))),cbind(diDF_s, dist = rep("Spearman",nrow(diDF_s))))
ggplot(diDF_s,aes(x=k,y=di,group= dist, color=dist,fill=group)) + geom_line() + geom_points()+ opts(title = "DI")

