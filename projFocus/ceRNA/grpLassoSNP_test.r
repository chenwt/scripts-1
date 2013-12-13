#!/usr/bin/Rscript
#Author: Jing He
#input: <file: expression and snp gt for one gene > 
#output: <file: gene: snp with weight 
#Usage: Rscript grplassoSNP.r input
#Description: require package irr, grpreg

#getting command line parameters
plotflag    = 0 # "0  for no plotting "
args        = commandArgs(TRUE)
if (is.null(args) || length(args) < 3){
  print("Please provide parameters in right order!")
  print("Usage: Rscript grplassoSNP_v1.r <inputsnp> <inputexp> <inputcnv> <1/0>")
  stop("Stop since no parameters")
}else{
  #print(args)
}

setwd(system("pwd",intern=T))
cwd         = getwd()
inputsnp    = args[1]
inputexp    = args[2]
inputcnv    = args[3]
plotflag    = args[4]
outputcoeff = paste(cwd,"/grplasso_coeff_",sep="")
# print(paste("current working directory",getwd()))
# print(paste("Current working directory:",system("pwd",intern=T)))

##-----------------------test-##----------------------------
#setwd("/Volumes/ac_lab/jh3283/scripts/projFocus/ceRNA/test/")
#setwd("/Volumes//ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test")
#setwd("/Users/jh3283/projFocus")
# inputsnp    = "input_test_reg_snp.mat"
# inputexp    = "input_test_reg_exp.mat"
# inputcnv    = "input_test_reg_cnv.mat"
# outputcoeff = "grplasso_coeff_"
##----------------------------##----------------------------

##-------useful functions
geneExpSigmoid  = function(x){
  #function to normalize and transform gene expression using sigmoid function
  x = sapply(x,as.numeric)
  y = (x - mean(x)) / sd(x)
  return(1/(1+exp(-y)))
}

##---svd
getSVD    = function(kcval){
  kcval_svdObj  = svd(kcval)
  temp          = cumsum(kcval_svdObj$d)
  cutoff_k      = min(which(temp > 0.8 * sum(kcval_svdObj$d),arr.ind=T))
  return(kcval_svdObj$u[,1:cutoff_k])
}


#--selecting optmial k using AIC 
bestK = function(kcval_svd){
  Ks          = 1: floor(nrow(kcval_svd)/3)
  aic         = (nrow(kcval_svd) - 1) * sum(apply(kcval_svd,2,var))
  aicdiff     = (nrow(kcval_svd) - 2) * sum(apply(kcval_svd,2,var))

  ###TODO tune lambda value for different data
  # cnt_lambda  = 0
  # bestKs      = NULL
  # for ( lambda in seq(0.5,1.5,by=0.02)){
  #   cnt_lambda = cnt_lambda + 1 
  #   for (i in Ks) {
  #       set.seed(12345)
  #       km      = kmeans(kcval_svd, centers=i)
  #       #AIC lambda = M (dimension)
  #       aic[i]  = km$tot.withinss + lambda * i
  #       bestKs  = c(bestKs, which.min(aic))
  #     }
  #   plot(aic,xlab=lambda)
  # }
  # plot(bestKs)
  # }
  lambda = 0.95
    for (i in Ks) {
        set.seed(12345)
        km      = kmeans(kcval_svd, centers=i)
        #AIC lambda = M (dimension)
        aic[i]  = km$tot.withinss + lambda * i
      }
  
  k_best  = which.min(aic)

  # # the first decreasing
  # for (i in 2:length(aic)) {
  #   aicdiff[i-1]   = (aic[i-1] - aic[i]) / aic[i-1]
  #   which
  #   if (aicdiff[i-1] > 0) {
  #     k_best  = i
  #     break
  #     }
  # }
  
  # Another method: the biggest r_square
  # r2=NULL
  # for (i in 3:(max(Ks)-3)) {
  # a       = lm(wss[1:i] ~ Ks[1:i])
  # b       = lm(wss[i+1:max(Ks)] ~ Ks[i+1:max(Ks)])
  # rsqd    = as.numeric(summary.lm(a)[8]) + as.numeric(summary.lm(b)[8])
  # r2      = c(r2, rsqd)}
  # k_best  = 3 + which(r2==max(r2))

  return(k_best=k_best)
}


getGroup = function(kcval_svd){
  k_best                = bestK(kcval_svd)
  colnames(kcval_svd)   = colnames(kcval)[1:ncol(kcval_svd)]
  row.names(kcval_svd)  = colnames(kcval)
  group                 = kmeans(x=kcval_svd,centers=k_best)$cluster
  return(unlist(group))
}
##----------------------------


# #---------loading snp data
# print("loading data")
# dataSnp             = read.table(inputsnp,header =T)
# rownames(dataSnp)   = dataSnp[,2]
# dataSnp             = dataSnp[,-c(1:2)]


# ##load cnv data
# dataCnv             = read.table(inputcnv,header =T)
# rownames(dataCnv)   = dataCnv[,2]
# dataCnv             = dataCnv[,-c(1:2)]

# ##load exp data
# dataExp             = read.table(inputexp,header =T)
# rownames(dataExp)   = dataExp[,1]
# genename            = dataExp[,1] 
# dataExp             = sapply(dataExp[,-c(1:4)],
#                             function(x){as.numeric(as.character(x))})


# ##------
# dataExpSg = geneExpSigmoid(dataExp)
# data_merge = t(rbind(dataExpSg,dataCnv,dataSnp))

# ##--------------------------------Grouping_of_variables
# library(irr)
# print("grouping of variables...")
# kcval = matrix(nrow=nrow(dataSnp), ncol= nrow(dataSnp))
# #kcp   = matrix(nrow=nrow(dataSnp), ncol= nrow(dataSnp))
# ###--TODO make it a parallel computing
# for ( i in 1:nrow(dataSnp))  
#   for ( j in i:nrow(dataSnp))    
#   {
#     kp2 = kappa2(ratings= t(rbind(dataSnp[i,],dataSnp[j,])) )
#     kcval[i,j]  = kcval[j,i] = kp2$value
#     #kcp[i,j]   = kcp[j,i]   = kp2$p.value 
#   }
# colnames(kcval) = row.names(dataSnp)
# rownames(kcval) = row.names(dataSnp)

# kcval_svd       = getSVD(kcval)
# group           = getGroup(kcval_svd)

# if (plotflag == 1){
#   pdf(paste(outputcoeff, genename,".pdf",sep=""))
#   plot(density(dataExp),main = "Gene Expression Density")
#   plot(density(dataExpSg), main = "Gene Expression Density after normalization")

# }

# ##########-----------model------------
# ###group lasso 
# #adding cnv variable
# #TODO: calculate cnv only residual, exclude confunding factor from CNV
# require(grpreg)
# group       = c(cnv = 1,group + 1)

# regfit      = function(X, y, group, plotflag=0){
#   #Desp:    function to perform regression 
#   #input:   datamatrix with first column as response, other as variables
#   #output:  coefficients beta, residual, r_square
#   #usage:   function internal called by regfitPermu
#   if (setequal(colnames(X), names(group))){
#         grpfit              = grpreg(X,y,group=group,penalty="grLasso")
#         grpfit_selected     = select(grpfit,'AIC')
#         if (plotflag == 1){
#           plot(grpfit)
#           plot(grpfit_selected$IC,xlab="Index of Lambda",ylab="IC")
#         }
#         beta                = grpfit_selected$beta[-1]
#         residual            = y - X %*% as.matrix(beta)
#         residual_sum_squre  = sum((residual)^2)
#         reg_sum_square      = sum((X %*% as.matrix(beta) - mean(y))^2)
#         total_sum_square    = reg_sum_square + residual_sum_squre
#         # r_square            = 1 - residual_sum_squre / total_sum_square #old one
#         r_square            = 1 - reg_sum_square / total_sum_square  # new one
#         return(list(beta    = beta,
#                   residual  = residual, 
#                   r_square  = r_square))
#   }else { print("ERROR! group not match data!")
#   }
# }


# regfitPermu       =   function(data_merge,group,nperm, plotflag=0) {
#   #Desp: function to calculate the significance of one fitting
#   #input: data_merge, group, and number of permuatation, plotflag 1/0
#   #output: p-value, plot of p-value significance depends on plotflag
#   y               =  data_merge[,1]  
#   X               =  as.matrix(data_merge[,-1])
#   fit_grpreg      = regfit(X, y, group, plotflag)

#   #permutation
#   count_perm      = 0 
#   r_square_perm   = 0.0
#   while (count_perm < as.numeric(nperm) ) {
#     y_permu             =   y + sample(fit_grpreg$residual,size=length(y),replace=TRUE)
#     y_permu             =   (y_permu - min(y_permu)) / (max(y_permu) - min(y_permu))
#     r_square_perm       =   c(r_square_perm, regfit(X, y_permu, group)$r_square) 
#     count_perm          =   count_perm  + 1
#   }
#   r_square_perm   = r_square_perm[-1]
#   pvalue          = length(r_square_perm[which(r_square_perm > fit_grpreg$r_square)]) / nperm

#   if (plotflag == 1){
#       plot(density(r_square_perm),
#            xlab = "",
#            sub  = paste("pvalue =",pvalue, "\n nperm =",nperm),
#            main = "Permutation results")
#       abline(v=fit_grpreg$r_square,col="red")
#   }
#     return(list(pvalue = pvalue,
#                   beta = fit_grpreg$beta))
  
# }

# print("Doing regression...")
# nperm     = 1000
# fitpermut = regfitPermu(data_merge,group,nperm, plotflag)
# # system.time(regfitPermu(data_merge,group,nperm=1000))
# # 80s for 1000 permutation

# ##output
# print("writing output...")
# outputcoeff = paste(outputcoeff,genename,".txt",sep="")
# write.table(t(as.matrix(c(paste("gene","npermu","pvalue",sep=":"),paste(genename, nperm, fitpermut$pvalue,sep=":")))),
#             outputcoeff,
#             quote=F,col.names=F,sep="\t",row.names = F)
# write.table(as.matrix(sort(fitpermut$beta)),
#             outputcoeff,
#             append=T,
#             col.names=F,quote=F,sep="\t")
