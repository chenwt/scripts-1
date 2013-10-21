myDunnIndex <- function(dataMA,k=3){
    mod <- kmeans(dataMA,k)
    cCenters <- mod$centers # matrix
    cSize <- mod$size #integer
    cDist <- matrix(NA,nrow = k,ncol = k)
    cMax <- list()

    for (i in 1:k){
        cMax[[i]] <- dataMA[names(which(mod$cluster == i)),] 
    }
    for ( i in 1:k){
        for( j in 1:k){
        cDist[i,j] <- myCdist(cMax[[i]],cMax[[j]])
        }
    }


    delta <- list()
    for (i in 1:k){
        delta[[i]] <- sum(myEdist(cMax[[i]],cCenters[i,])) / cSize[i]
    }

    delta
    di <- matrix(NA, nrow=k,ncol=k)
    for( i in 1:k){
        for ( j in 1:k)
        {
            if( delta[[i]] > 0 || delta[[j]] > 0 ){
                ditmp<- myCdist(cMax[[i]],cMax[[j]])/max(delta[[i]],delta[[j]])
                di[i,j] <- ditmp
            }else {
                di[i,j] <- 10000 
            }
        }
    }
    diRe <- min(di)
    return(diRe)
}


myCdist <- function(x,y){

    nr <- nrow(x)
    nc <- nrow(y)
    if(is.null(nc)) nc <- 1
    if(is.null(nr)) nr <- 1
    dmax <- 0
    if(nc > 1 && nr > 1){
        for ( i in (1:nc)){
          dmax <- max(max(sqrt(rowSums((x-y[i,])^2))),dmax)
        }
    } else {
        if( nc/nr > 1) {tmp <- x;x <- y; y <- tmp}
        dmax <- max(myEdist(x,y))
    }
    return(dmax)
}

myEdist <- function(x,y){
    nr <- nrow(x)
    if(is.null(nr)) nr <- 1
    # for ( i in (1:nc)){
    if(nr > 1) {
        d <- matrix(NA,nrow=nr,ncol=1)
        d[,i] <- sqrt(rowSums((x-y)^2))
    }else{
         d <- sqrt(sum((x-y)^2))
    }
    return(d)
}



myDiBatch <- function(dataMA,mink,maxk){
    diDF <- data.frame(cbind(k=mink:maxk,di=NA))

    for(i in mink:maxk)
     diDF[(i-mink+1),2]<- as.numeric(myDunnIndex(dataMA，i))
    return(diDF)
}

diDF <- myDiBatch(exp177MA,2,60)
ggplot(diDF,aes(x=k,y=di)) + geom_line(color="blue")


ggplot(diDF) + geom_line(aes(x=k, y=value, colour=variable)) +
  scale_colour_manual(values="blue")


