#!/bin/Rscirpts
#Author: Jing He
#Date: Apr.2nd, 2013
#Last Update: Apr.25th, 2013

require("mixtools")

#attach(faithful)
if(Sys.info()["sysname"] == "Darwin") {
	setwd("/Volumes/ys_lab_scratch/jh3283/school/compGenomic/mixtools/")
	print("work on IMAC!")
}else {
	setwd("~/SCRATCH/school/compGenomic/mixtools")
	print("Work On TITAN!")
}


getData <- function(filename) {
	data <- read.table(filename,header=F)
	# ac <- data[,3:4]
	# row.names(ac) <- unlist(paste(data[,1],data[,2],sep="_"))
	data$BAF <- data[,3]/data[,4]
	colnames(data) <- c("Chr","Pos","Alt","Ref","BAF")
	return(data)
}

##------------------------------for one sample 
myGaussianEMPlot <- function(data,out,lambda,mu){
	pdf(paste(out,"_2GaussianEM.pdf",sep=""))
	baf1 <- normalmixEM(data,lambda=lambda,mu=mu)
	plot(density(data),type="l")
	plot(baf1, density = TRUE, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8,
				main2 = paste(out,"BAF of genom sites",sep=""), xlab2 = "BAF")
	dev.off()
}
myGaussianEM <- function(data,lambda,mu){
	baf1 <- normalmixEM(data,lambda=lambda,mu=mu)
	params <- data.frame(baf1[c("lambda", "mu", "sigma")])
	return(params)
}

##------------------------------for multiple sample
ac53 <- getData("ac53.freq")
data <- ac53$BAF

fns <- dir("tumorFreq")[grep("freq", dir())]
ac53<- getData(fns[1])
ac54<- getData(fns[2])
ac55<- getData(fns[3])
ac56<- getData(fns[4])
acAll <- merge (merge(ac53,ac54,by=c("Chr","Pos")),ac55,by=c("Chr","Pos"))
acAll <- acAll[,c("BAF.x","BAF.y","BAF")]
plot(density(acAll[,1]),type="l",main="")
lines(density(acAll[,2]),col="blue")
lines(density(acAll[,3]),col="green")
title("Three Tumor Sample density Plot")

p1 <- myGaussianEM(ac53$BAF,"ac53")
p2 <- myGaussianEM(ac54$BAF,"ac54")
p3 <- myGaussianEM(ac55$BAF,"ac55")
p4 <- myGaussianEM(ac56$BAF,"ac56")

cutpts<-c(0.3, 0.5, 0.8)
mulData <- makemultdata(acAll, cuts = cutpts)
set.seed(15)
theta4 <- matrix(runif(nrow(acAll)*3), ncol = 3)
theta3 <- theta4[1:3,]
mult3 <- multmixEM(mulData, lambda=c(0.1,0.4,0.5) )


##------------------------------updating 
fns <- dir("tumorFreq")
tumor <- list()
for (i in 1:length(fns)) tumor[[i]] <- getData(paste("tumorFreq/",fns[i],sep=""))

resParaLS <- list()
for (i in 1: length(tumor)) 
 {resParaLS[[i]] <- myGaussianEM(tumor[[i]]$BAF,lambda=.5,mu=c(.3,.8))}

resMX <- cbind(unlist(lapply(resParaLS, function(x){ x$lambda[which.min(x$mu)]}))
	,unlist(lapply(resParaLS, function(x){ x$lambda[which.max(x$mu)]})))

heatmap(resMX
, scale=c("row", "column", "none")
)
match.by <- c("Chr","Pos")
tumorMrg <- Reduce(function(...) merge(..., by=match.by), tumor)
colnames(tumorMrg) <- c("Chr","Pos",paste(rep(gsub(".freq","",fns),each=3),rep(c("Alt","Ref","BAF"),length(tumor)),sep="."))

##------------------------------normal data 
# normal <- getData("../vcfs/ac1.freq")
normal <- getData("ac1.freq")
cosmic <- read.table("cosmic.pos",header=F)
colnames(cosmic) <- c("Chr","Pos")
# nor <- merge(normal,cosmic,by=c("Chr","Pos")) less 

myGaussianEM(nor$Alt/(nor$Ref + nor$Alt))
myGaussianEMPlot(nor$Alt/(nor$Ref + nor$Alt),"ac1",lambda =.5, mu=c(.1,1))
normal$BAF <- normal$Alt/(normal$Ref + normal$Alt)
baf2 <- normalmixEM(normal$Alt/(normal$Ref + normal$Alt),lambda=.5,mu=c(0.1,.5))

normalmixEM(no$BAF,k=2,maxit=100,epsilon=0.01)

plot(hist(normal$BAF[normal$BAF>0],breaks=200),col="grey",border="grey",freq=FALSE,
xlab="",main="")
lines(density(normal$BAF[normal$BAF>0]),lty=2)



plot.gaussion.components <- function(mixture,component.number,...) {
	curve(mixture$lambda[component.number] *
	dnorm(x,mean=mixture$mu[component.number],
	sd=mixture$sigma[component.number]), add=TRUE, ...)
}

plot(hist(tumor[[1]]$BAF,breaks=200),col="grey",border="grey",freq=FALSE,
xlab="",main="")
lines(density(tumor[[1]]$BAF),lty=2)
tm1 <- normalmixEM(tumor[[1]]$BAF,lambda=.5,mu=c(.3,.8))
sapply(1:2,plot.gaussion.components,mixture=normalmixEM(tm1))

plot(density(tumor[[1]]$BAF),col="red",type="l", main="", xlab="",ylab="")
lines(density(tumor[[2]]$BAF),col="red",pch="*")
lines(density(tumor[[3]]$BAF),col="red",pch="*")
lines(density(tumor[[4]]$BAF),col="red",pch="*")
lines(density(nor$Alt/(nor$Ref + nor$Alt)),col="green")
title(main="density of BAF for different samples", xlab="BAF",y="density")

##------------------------------ cross validation to select K


dnormalmix <- function(x,mixture,log=FALSE) {
	lambda <- mixture$lambda
	k <- length(lambda)
# Calculate share of likelihood for all data for one component
	like.component <- function(x,component) {
		lambda[component]*dnorm(x,mean=mixture$mu[component],
		sd=mixture$sigma[component])
	}
# Create array with likelihood shares from all components over all data
	likes <- sapply(1:k,like.component,x=x)
# Add up contributions from components
	d <- rowSums(likes)
	if (log) {
		d <- log(d)
	}
	return(d)
}
loglike.normalmix <- function(x,mixture) {
loglike <- dnormalmix(x,mixture,log=TRUE)
return(sum(loglike))
}

# check 
loglike.normalmix(tumor[[1]]$BAF,mixture=tm1)


t1data <- tumor[[1]]$BAF
n <- length(t1data)
data.points <- 1:n
data.points <- sample(data.points) # Permute randomly
train <- data.points[1:floor(9 * n /10)] # First random half is training
test <- data.points[-(1:floor(n/10))]
candidate.component.numbers <- 2:10
loglikes <- vector(length=1+length(candidate.component.numbers))
mu <- mean(t1data[train])
sigma <- sd(t1data[train])*sqrt((n-1)/n) # MLE of standard deviation
loglikes[1] <- sum(dnorm(t1data[test],mu,sigma,log=TRUE))

for (k in candidate.component.numbers) {
	mixture <- normalmixEM(t1data[train],k=k,maxit=500,epsilon=1e-2)
	loglikes[k] <- loglike.normalmix(t1data[test],mixture=mixture)
}

loglikes

plot(x=1:10, y=loglikes,xlab="Number of mixture components",
ylab="Log-likelihood on testing data", main="t1")
points(x=4,y=loglikes[[4]],col="red",pch=19)

## plot components
t1data.k4 <- normalmixEM(t1data,k=9,maxit=500,epsilon=1e-2)
plot(hist(t1data,breaks=101),col="grey",border="grey",freq=FALSE,
xlab="",main="")
lines(density(t1data),lty=2)
# sapply(1:4,plot.gaussion.components,mixture=t1data.k9,col="red")
mycol <- c("blue","red","green","yellow")
for (i in 1:4) {
	plot.gaussion.components(t1data.k9,component.number=i,col=mycol[i])
}



##------------------------------calibration 

pnormmix <- function(x,mixture) {
	lambda <- mixture$lambda
	k <- length(lambda)
	pnorm.from.mix <- function(x,component) {
	lambda[component]*pnorm(x,mean=mixture$mu[component],
	sd=mixture$sigma[component])
	}
	pnorms <- sapply(1:k,pnorm.from.mix,x=x)
	return(rowSums(pnorms))
}

t1data.k2 <- normalmixEM(t1data,k=2,maxit=500,epsilon=1e-2)
par(mfrow=c(2,1))
distinct.t1data <- sort(unique(t1data))
tcdfs <- pnormmix(distinct.t1data,mixture=t1data.k2)
ecdfs <- ecdf(t1data)(distinct.t1data)
plot(tcdfs,ecdfs,pch = 20, cex=.3
	, xlab="Theoretical CDF",ylab="Empirical CDF"
	,xlim=c(0,0.2),ylim=c(0,0.2)
	,main="two components",col="blue")
abline(0,1)

distinct.t1data <- sort(unique(t1data))
tcdfs <- pnormmix(distinct.t1data,mixture=t1data.k4)
ecdfs <- ecdf(t1data)(distinct.t1data)
plot(tcdfs,ecdfs,pch = 20, cex=.3,
	xlab="Theoretical CDF",ylab="Empirical CDF"
	,xlim=c(0,0.2),ylim=c(0,.2)
	,main="four components",col="blue")
abline(0,1)



##------------------------------visulize the data
plot(rbind(p1$lambda,p2$lambda,p3$lambda)[,1],rbind(p1$lambda,p2$lambda,p3$lambda)[,2]
		,xlim=c(0,1),ylim=c(0,1))







