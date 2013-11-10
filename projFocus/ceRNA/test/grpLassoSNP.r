#!/usr/bin/Rscript
#Author: Jing He
#input: <file: expression and snp gt for one gene > 
#output: <file: gene: snp with weight 
#Usage: Rscript grplassoSNP.r input
#Description: 


setwd("/Volumes/ac_lab/jh3283/scripts/projFocus/ceRNA/test/")

## Linear regression
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)

## Logistic regression
y <- birthwt.grpreg$low
fit <- grpreg(X,y,group,penalty="grLasso", family="binomial")
plot(fit)
fit <- grpreg(X,y,group,penalty="grMCP", family="binomial")
plot(fit)
fit <- grpreg(X,y,group,penalty="grSCAD", family="binomial")
plot(fit)
fit <- grpreg(X,y,group,penalty="gMCP", family="binomial")
plot(fit)
select(fit,"BIC")
