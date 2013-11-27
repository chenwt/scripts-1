#!/usr/bin/Rscript
#Author: Jing He
#input: <file: expression and snp gt for one gene > 
#output: <file: gene: snp with weight 
#Usage: Rscript grplassoSNP.r input
#Description: 

# args = commandArgs(TRUE)
#  if (is.null(args)){
#       print("Please provide parameters")
#    exit
#     }else{
#          print(args)
#     }


#setwd("/Volumes/ac_lab/jh3283/scripts/projFocus/ceRNA/test/")
setwd("/Volumes//ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test")
##----------------------------
#setting parameters
# inputmut = args[1]
# inputmut = "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test_knowBRCA.gene_snp_meth_cnv"
# inputexp = "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test_knowBRCA.gene_exp" 

inputsnp = "input_test_reg_snp.mat"
inputexp = "input_test_reg_exp.mat"
inputcnv = "input_test_reg_cnv.mat"

##----------------------------
#loading snp data
dataSnp            = read.table(inputsnp,header =T)
rownames(dataSnp)  = dataSnp[,2]
dataSnp = dataSnp[,-c(1:2)]


##load cnv data
dataCnv = read.table(inputcnv,header =T)
rownames(dataCnv)  = dataCnv[,2]
dataCnv = dataCnv[,-c(1:2)]

##load exp data
dataExp           = read.table(inputexp,header =T)
rownames(dataExp) = dataExp[,1]
genename  = dataExp[,1] 
dataExp   = sapply(dataExp[,-c(1:4)],function(x){as.numeric(as.character(x))})
plot(density(dataExp))


##-------useful functions
geneExpSigmoid = function(x){
  #function to normalize and transform gene expression using sigmoid function
  x = sapply(x,as.numeric)
  y = (x - mean(x)) / sd(x)
  return(1/(1+exp(-y)))
}

##------
dataExpSg = geneExpSigmoid(dataExp)
plot(density(dataExpSg))
data_merge = t(rbind(dataExpSg,dataCnv,dataSnp))



##------------Grouping_of_variables
library(irr)
dataSnp_new <- data.frame(R1="",R2="",stringsAsFactors=FALSE)
for(x in 1:ncol(dataSnp)) {
    dataSnp_new <- rbind(dataSnp_new,rep(rownames(dataSnp),dataSnp[,x]))
}
rm(x)
dataSnp_new <- dataSnp_new[-1,]
kappa2(ratings=dataSnp_new )
kcval = matrix(nrow=nrow(dataSnp), ncol= nrow(dataSnp))
kcp   = matrix(nrow=nrow(dataSnp), ncol= nrow(dataSnp))
###----TODO make it a parallel computing
for ( i in 1:nrow(dataSnp))  
  for ( j in i:nrow(dataSnp))    
  {
     kp2 = kappa2(ratings= t(rbind(dataSnp[i,],dataSnp[j,])) )
     kcval[i,j]  = kcval[j,i] = kp2$value
     #kcp[i,j]   = kcp[j,i]   = kp2$p.value 
  }

kcval_svd = svd(kcval)
save.image("grpLassoSNP.RData")
Ks = 1: floor(nrow(kcval) / 2)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in Ks) wss[i] <- sum(kmeans(mydata, 
     centers=i)$withinss)

### function to do logistic regression to one gene data

require("glmnet")

gla_cv  =     cv.glmnet(y = data_merge[,1], x = data_merge[,-1], 
                nfolds = 10,alpha = 1,type.gaussian="covariance" ) 

for (alpha in seq(0,1,by=0.05)) {
    gla_fit <- glmnet(y = data_merge[,1], x = data_merge[,-1], 
                  alpha = alpha,
                  type.multinomial = "grouped")
    plot(gla_fit,
         main = paste("alpha = ",alpha) )
}
plot(gla_fit)

gla$lambda.1se
gla_best  =   gla$glmnet.fit
gla.coef  =   coef(gla$glmnet.fit, s = gla$lambda.1se)

lambda    =   gla_cv$lambda.min
gla_fit   =   glmnet(y = data_merge[,1], x = data_merge[,-1], 
                     alpha = alpha,
                     type.gaussian="covariance" ,lambda=lambda)
coeff = as.matrix(t(coef(gla_fit)))
coeff = coeff[order(coeff)]
plot(x=1:length(coeff), y = coeff,
     col='blue',
     main = paste("lambda =",round(lambda,5)) )
text(x=1:length(coeff), y = coeff,labels=colnames(coeff))
predict(gla_fit,"coefficients")
plot(gla_fit)
gla.coef[which(gla.coef != 0)]

#  install.packages("lars")
require(lars)
reg_lasso    =   lars(y = data_merge[,1], x = data_merge[,-1],
                   type = "lasso",intercept = TRUE)
reg_lasso_cv =   cv.lars(y = data_merge[,1], x = data_merge[,-1],
                     type = "lasso",trace = TRUE)
#plot(reg_lasso,trace = F)

# require("glmnetcr")
# glmnet.fit <- glmnet.cr(y = data_merge[,1], x = data_merge[,-1],alpha = 1)
# AIC <- select.glmnet.cr(glmnet.fit, which="AIC")
# fitted(glmnet.fit, s=AIC)


# ######----------------------------
# #this part under development

# expSCReg = function (mdata,edata){
#   mdmat = t(mdata) 
#   dmat = data.frame(edata,mdmat)
#   glm_out = glm(edata ~ mdmat , family=binomial(logit)) 
#   lm_logit = lm(edata ~ mdmat)
#   pdf("lasso_logit.pdf")

#   plot(glm_out)
#   plot(lm_linear)
#   plot(lm_logit)



#   reg_lasso_cv = cv.lars(x=mdmat, y = edata, K = 10, plot.it = TRUE)
#   # reg_lasso = lars(x = mdmat, y = t(edata) ,type = "lasso")
#   reg_lasso_cv = cv.lars(x=mdmat, y = t(edata), K = 10, plot.it = TRUE)
#   best <- reg_lasso_cv$index[which.min(reg_lasso_cv$cv)]
#   coef <- coef.lars(reg_lasso, mode = "fraction", s = best) 
#   coef[coef != 0] 




#   # plot(reg_lasso,breaks=F)
  
#   plot(gla.best)
#   plot(reg_lasso,breaks =F)
#   boxplot(reg_lasso$beta,xlab="genomic_mutations",y="coefficients",pch = "*")

#   dev.off()

# } 
# ##----------------------------development-end--------------
