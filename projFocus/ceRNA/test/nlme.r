##!/usr/bin/Rscript
#Author: Jing He
#Date: 
#Last Updated:
#Usage:
#Description: 

require("lme4")
require("nlme")

##----------------------------
#test

lmm.data <- read.table("http://www.unt.edu/rss/class/Jon/R_SC/Module9/lmm.data.txt", 
	header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE) 
lmm.2 <- lmer(formula = extro ~ open + agree + social + class + (1 |school/class), data=lmm.data, family = gaussian, REML = TRUE, verbose = FALSE) 

##----------------------------
#prepare data 


##----------------------------
#setting parameter


##----------------------------
#fitting model


##----------------------------
#output


##----------------------------
#clean up

> Th . start <- c( lKe = -2.5 , lKa = 0.5 , lCl = -3)
> nm1 <- nlmer ( conc ~ SSfol ( Dose , Time , lKe , lKa , lCl ) ~
+ 0+ lKe + lKa + lCl +(0+ lKe | Subject )+(0+ lKa | Subject )
+ +(0+ lCl | Subject ), nAGQ =0 , Theoph ,
+ start = Th . start , verbose = TRUE )


fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
            data = Loblolly,
            fixed = Asym + R0 + lrc ~ 1,
            random = Asym ~ 1,
            start = c(Asym = 103, R0 = -8.5, lrc = -3.3))
summary(fm1)
fm2 <- update(fm1, random = pdDiag(Asym + lrc ~ 1))
summary(fm2)