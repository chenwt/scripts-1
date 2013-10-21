#!/bin/Rscirpts
#Author: Jing He
#Date: Apr.2nd, 2013
#Last Update: Apr.7th, 2013

require("mixtools")

attach(faithful)
if(Sys.info()["sysname"] == "Darwin") {
	setwd("/Volumes/ys_lab_scratch/jh3283/school/compGenomic/mixtools/")
	print("work on IMAC!")
}else {
	setwd("~/SCRATCH/school/compGenomic/mixtools")
	print("Work On TITAN!")
}


getData <- function(filename){
	data <- read.table(filename,header=F)
	data$BAF <- as.numeric(data[,4])/(as.numeric(data[,3]) + as.numeric(data[,4]))
	colnames(data) <- c("Chr","Pos","Ref","Alt","BAF")
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

# ##------------------------------for multiple sample

plot.gaussion.components <- function(mixture,component.number,...) {
	curve(mixture$lambda[component.number] *
	dnorm(x,mean=mixture$mu[component.number],
	sd=mixture$sigma[component.number]), add=TRUE, ...)
}

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

##------------------------------Maximum likelihood for different candidate component numbers
myMultiMixMLE <- function(baf,kmix){
	n <- length(baf)
	dpts <- 1:n
	dpts <- sample(dpts) # Permute randomly
	train <- dpts[1:floor(9 * n /10)] # First random half is training
	test <- dpts[-(1:floor(n/10))]

	numMixs <- 2:kmix
	loglikes <- vector(length=1+length(numMixs))
	mu <- mean(baf[train])
	sigma <- sd(baf[train]) * sqrt((n - 1)/n) # MLE of standard deviation
	loglikes[1] <- sum(dnorm(baf[test],mu,sigma,log=TRUE))

	for (k in numMixs) {
		mixture <- normalmixEM(baf[train],k=k,maxit=500,epsilon=1e-2)
		loglikes[k] <- loglike.normalmix(baf[test],mixture=mixture)
	}
	loglikes
}

getKMix <- function(baf,loglikes) {
	lagdiff <- diff(loglikes,lag=1)
	pll <- loglikes[which.min(lagdiff)]
	km <- which.min(lagdiff)

	pdf(paste(fnames,"mle.pdf",sep=""))
	plot(x=1:kmix, y=loglikes,xlab="Number of mixture components",
		ylab="Log-likelihood on testing data", main=fnames)
	points(x=km,y=pll,col="red",pch=19)
	## plot components
	KmixFit <- normalmixEM(baf,k=km,maxit=500,epsilon=1e-2)
	plot(hist(baf,breaks=101),col="grey",
		border="grey",freq=FALSE,
		xlab="",main="")
	lines(density(baf),lty=2)
	# sapply(1:4,plot.gaussion.components,mixture=KmixFit,col="red")
	# mycol <- c("blue","red","green","yellow")
	for (i in 1:km) {
		plot.gaussion.components(KmixFit,component.number=i,col=2*km[i])
	}
	dev.off()
	
	return(KmixFit)
}

