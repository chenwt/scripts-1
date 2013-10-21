#!/bin/Rscirpts
#Author: Jing He
#Date: Apr.2, 2013
#Last Update: Apr. 7, 2013
#parameters: 
#example: Rscript 

###------------------------------header------------------------------
# args <- commandArgs(TRUE)

# if (is.null(args)){
# 	print("Please provide parameters")
# 	exit
# }else{
# 	print(args)
# }

# ###------------------------------start coding##------------------------------

# require("DPpackage")


# library("coda")
# coda.obj <- mcmc.list(
# chain1 = mcmc(fit1$save.state$thetasave[,1]),
# chain2 = mcmc(fit2$save.state$thetasave[,1]),
# chain3 = mcmc(fit3$save.state$thetasave[,1]))
# gelman.diag(coda.obj, transform = TRUE)




# dtrue <- function(grid, x) {
# 	exp(-2 * x) * dnorm(grid, mean = x, sd = sqrt(0.01)) +
# 	(1 - exp(-2 * x)) * dnorm(grid, mean = x^4, sd = sqrt(0.04))
# }
# mtrue <- function(x) exp(-2 * x) * x + (1 - exp(-2 * x)) * x^4
# 	# The data were simulated using the following code:
# set.seed(0)
# nrec <- 500
# x <- runif(nrec)
# y1 <- x + rnorm(nrec, 0, sqrt(0.01))
# y2 <- x^4 + rnorm(nrec, 0, sqrt(0.04))
# u <- runif(nrec)
# prob <- exp(-2 * x)
# y <- ifelse(u < prob, y1, y2)

#  w <- cbind(y, x)
#  wbar <- apply(w, 2, mean)
#  wcov <- var(w)
#  prior <- list(a0 = 10, b0 = 1, nu1 = 4, nu2 = 4, s2 = 0.5 * wcov,
# m2 = wbar, psiinv2 = 2 * solve(wcov), tau1 = 6.01, tau2 = 3.01)

# mcmc <- list(nburn = 5000, nsave = 5000, nskip = 3, ndisplay = 1000)

# fitWDDP <- DPcdensity(y = y, x = x, xpred = seq(0, 1, 0.02),
# 	ngrid = 100, compute.band = TRUE, type.band = "HPD",
# 	prior = prior, mcmc = mcmc, state = NULL, status = TRUE)




## Not run: 
    # Respiratory Data Example

      data(indon)
      attach(indon)

      baseage2<-baseage**2
      follow<-age-baseage
      follow2<-follow**2 

    # Prior information

      prior<-list(alpha=1,
                  nu0=4.01,
                  tinv=diag(1,1),
                  nub=4.01,
                  tbinv=diag(1,1),
                  mb=rep(0,1),
                  Sb=diag(1000,1),
                  beta0=rep(0,9),
                  Sbeta0=diag(1000,9))

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-5000
      nsave<-25000
      nskip<-20
      ndisplay<-1000
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Fit the Probit model
      fit1<-DPMglmm(fixed=infect~gender+height+cosv+sinv+xero+baseage+
                    baseage2+follow+follow2,
                    random=~1|id,family=binomial(probit),
                    prior=prior,mcmc=mcmc,state=state,status=TRUE)

    # Fit the Logit model
      fit2<-DPMglmm(fixed=infect~gender+height+cosv+sinv+xero+baseage+
                    baseage2+follow+follow2,random=~1|id,
                    family=binomial(logit),
                    prior=prior,mcmc=mcmc,state=state,status=TRUE)
      save.image("dppkg_test.Rdata")
    # Summary with HPD and Credibility intervals
      print(summary(fit1))
       print(summary(fit1,hpd=FALSE))

      summary(fit2)
      summary(fit2,hpd=FALSE)


    # Plot model parameters (to see the plots gradually set ask=TRUE)
      pdf("dp_test.pdf")
      plot(fit1,ask=FALSE)
      plot(fit1,ask=FALSE,nfigr=2,nfigc=2)	

    # Plot an specific model parameter (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="baseage")	
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="ncluster")	

      dev.off()
## End(Not run)

####################################
    # Univariate example
    ####################################

    # Data
      data(galaxy)
      galaxy <- data.frame(galaxy,speeds=galaxy$speed/1000) 
      attach(galaxy)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 1000
      nsave <- 10000
      nskip <- 10
      ndisplay <- 100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Example of Prior information 1
    # Fixing alpha, m1, and Psi1

      prior1 <- list(alpha=1,m1=rep(0,1),psiinv1=diag(0.5,1),nu1=4,
                     tau1=1,tau2=100)


    # Example of Prior information 2
    # Fixing alpha and m1

      prior2 <- list(alpha=1,m1=rep(0,1),psiinv2=solve(diag(0.5,1)),
                     nu1=4,nu2=4,tau1=1,tau2=100)


    # Example of Prior information 3
    # Fixing only alpha

      prior3 <- list(alpha=1,m2=rep(0,1),s2=diag(100000,1),
                   psiinv2=solve(diag(0.5,1)),
                   nu1=4,nu2=4,tau1=1,tau2=100)


    # Example of Prior information 4
    # Everything is random

      prior4 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1),
                   psiinv2=solve(diag(0.5,1)),
                   nu1=4,nu2=4,tau1=1,tau2=100)

    # Fit the models

      fit1.1 <- DPdensity(y=speeds,prior=prior1,mcmc=mcmc,
                          state=state,status=TRUE)
      fit1.2 <- DPdensity(y=speeds,prior=prior2,mcmc=mcmc,
                          state=state,status=TRUE)
      fit1.3 <- DPdensity(y=speeds,prior=prior3,mcmc=mcmc,
                          state=state,status=TRUE)
      fit1.4 <- DPdensity(y=speeds,prior=prior4,mcmc=mcmc,
                          state=state,status=TRUE)

    # Posterior means
      fit1.1
      fit1.2
      fit1.3
      fit1.4

    # Plot the estimated density
      plot(fit1.1,ask=FALSE)
      plot(fit1.2,ask=FALSE)
      plot(fit1.3,ask=FALSE)
      plot(fit1.4,ask=FALSE)

    # Extracting the density estimate
      cbind(fit1.1$x1,fit1.1$dens)
      cbind(fit1.2$x1,fit1.2$dens)
      cbind(fit1.3$x1,fit1.3$dens)
      cbind(fit1.4$x1,fit1.4$dens)
      
    # Plot the parameters (only prior 2 for illustration)
    # (to see the plots gradually set ask=TRUE)
      plot(fit1.2,ask=FALSE,output="param")

    # Plot the a specific parameters 
    # (to see the plots gradually set ask=TRUE)
      plot(fit1.2,ask=FALSE,output="param",param="psi1-speeds",
           nfigr=1,nfigc=2)

    # Extracting the posterior mean of the specific 
    # means and covariance matrices 
    # (only prior 2 for illustration)
      DPrandom(fit1.2) 

    # Ploting predictive information about the specific 
    # means and covariance matrices 
    # with HPD and Credibility intervals
    # (only prior 2 for illustration)
    # (to see the plots gradually set ask=TRUE)
      plot(DPrandom(fit1.2,predictive=TRUE),ask=FALSE)
      plot(DPrandom(fit1.2,predictive=TRUE),ask=FALSE,hpd=FALSE)

    # Ploting information about all the specific means 
    # and covariance matrices 
    # with HPD and Credibility intervals
    # (only prior 2 for illustration)
    # (to see the plots gradually set ask=TRUE)
      plot(DPrandom(fit1.2),ask=FALSE,hpd=FALSE)


    ####################################
    # Bivariate example
    ####################################

    # Data
      data(airquality)
      attach(airquality)

      ozone <- Ozone**(1/3)
      radiation <- Solar.R

    # Prior information

      s2 <- matrix(c(10000,0,0,1),ncol=2)
      m2 <- c(180,3)
      psiinv2 <- solve(matrix(c(10000,0,0,1),ncol=2))
     
      prior <- list(a0=1,b0=1/5,nu1=4,nu2=4,s2=s2,
                    m2=m2,psiinv2=psiinv2,tau1=0.01,tau2=0.01)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5000
      nsave <- 10000
      nskip <- 10
      ndisplay <- 1000
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Fit the model
      fit1 <- DPdensity(y=cbind(radiation,ozone),prior=prior,mcmc=mcmc,
                        state=state,status=TRUE,na.action=na.omit)

    # Plot the estimated density
      plot(fit1)

    # Extracting the density estimate
      fit1$x1
      fit1$x2
      fit1$dens

## End(Not run)





