#!/bin/Rscirpts
#Author: Jing He
#Date: Apr.2, 2013
#Last Update: 
#parameters: 
#example: Rscript 

###------------------------------header------------------------------
args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters")
	exit
}else{
	print(args)
}

###------------------------------start coding##------------------------------

require("DPpackage")
library("coda")
coda.obj <- mcmc.list(
chain1 = mcmc(fit1$save.state$thetasave[,1]),
chain2 = mcmc(fit2$save.state$thetasave[,1]),
chain3 = mcmc(fit3$save.state$thetasave[,1]))
gelman.diag(coda.obj, transform = TRUE)




dtrue <- function(grid, x) {
	exp(-2 * x) * dnorm(grid, mean = x, sd = sqrt(0.01)) +
	(1 - exp(-2 * x)) * dnorm(grid, mean = x^4, sd = sqrt(0.04))
}
mtrue <- function(x) exp(-2 * x) * x + (1 - exp(-2 * x)) * x^4
	# The data were simulated using the following code:
set.seed(0)
nrec <- 500
x <- runif(nrec)
y1 <- x + rnorm(nrec, 0, sqrt(0.01))
y2 <- x^4 + rnorm(nrec, 0, sqrt(0.04))
u <- runif(nrec)
prob <- exp(-2 * x)
y <- ifelse(u < prob, y1, y2)

 w <- cbind(y, x)
 wbar <- apply(w, 2, mean)
 wcov <- var(w)
 prior <- list(a0 = 10, b0 = 1, nu1 = 4, nu2 = 4, s2 = 0.5 * wcov,
m2 = wbar, psiinv2 = 2 * solve(wcov), tau1 = 6.01, tau2 = 3.01)

mcmc <- list(nburn = 5000, nsave = 5000, nskip = 3, ndisplay = 1000)

fitWDDP <- DPcdensity(y = y, x = x, xpred = seq(0, 1, 0.02),
	ngrid = 100, compute.band = TRUE, type.band = "HPD",
	prior = prior, mcmc = mcmc, state = NULL, status = TRUE)





