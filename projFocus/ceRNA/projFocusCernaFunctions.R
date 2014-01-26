##-------useful functions-------------------
getArgs = function(){
  #base function used to get commandline argument specified by -- or by sequence
  args = commandArgs(trailingOnly = TRUE)
  hh <- paste(unlist(args),collapse=' ')
  if(length(grep("--",hh)) > 0){
  listOptions <- unlist(strsplit(hh,'--'))[-1]
  optionsArgs <- sapply(listOptions,function(x){
           unlist(strsplit(x, ' '))[-1]
          })
  optionsNames <- sapply(listOptions,function(x){
    option <-  unlist(strsplit(x, ' '))[1]
  })
  names(optionsArgs) <- unlist(optionsNames)
 }else{
  optionsArgs = args
 }
  return(optionsArgs)
}


getData = function(inputexp,inputsnp,inputcnv){
  ##load exp data
  dataExp             = read.table(inputexp,header =T)
  rownames(dataExp)   = dataExp[,1]
  genename            = dataExp[,1]
  print(class(dataExp))
  dataExp             = sapply(dataExp[,-c(1:4)],
                              function(x){as.numeric(as.character(x))})
  dataExpNorm           = normalize(dataExp)
  #---------loading snp data
  dataSnp             = read.table(inputsnp,header =T)
  rownames(dataSnp)   = dataSnp[,2]
  dataSnp             = dataSnp[,-c(1:2)]
  #snp data QC:1. delete all 0 snp
  dataSnp             = dataSnp[rowSums(dataSnp)>0,]

  if (nrow(dataSnp) > 0 ){
      ##load cnv data
      dataCnv             = read.table(inputcnv,header =T)
      rownames(dataCnv)   = dataCnv[,2]
      dataCnv             = dataCnv[,-c(1:2)]

      ##------
      data_merge  = t(rbind(dataExpNorm,dataCnv,dataSnp))
    }
  return(list(dataSnp=dataSnp,
              dataExp=dataExp,
              dataCnv=dataCnv,
              data_merge=data_merge,
              genename = genename))
}

getAllData = function(inputexp,inputsnp,inputcnv,inputsom){
  ##when all data avaiable
  #-------load exp data
  dataExp             = read.table(inputexp,header =T)
  rownames(dataExp)   = as.character(dataExp[,1])
  genename            = dataExp[,1] 
  dataExp             = sapply(dataExp[,-c(1:4)],
                              function(x){as.numeric(as.character(x))})
  dataExpNorm           = as.data.frame(t(normalize(dataExp)))

  #--------load snp data
  dataSnp             = read.table(inputsnp,header =T)
  rownames(dataSnp)   = dataSnp[,2]
  dataSnp             = dataSnp[,-c(1:2)]
  #snp data QC:1. delete all 0 snp
  dataSnp             = dataSnp[rowSums(dataSnp)>0,]

  if (nrow(dataSnp) > 0 ){
      ##load cnv data
      dataCnv             = read.table(inputcnv,header = T)
      rownames(dataCnv)   = dataCnv[,2]
      dataCnv             = dataCnv[,-c(1:2)]
      ##load somatic mutation data
      dataSom             = read.table(inputsom,header = T)
      rownames(dataSom)   = "som"
      dataSom             = dataSom[,-c(1:2)]
      ##------
      names(dataSom)      = sapply(names(dataSom),subStr1To19)
      names(dataCnv)      = sapply(names(dataCnv),subStr1To19)
      names(dataSnp)      = sapply(names(dataSnp),subStr1To19)
      names(dataExpNorm)  = sapply(names(dataExpNorm),subStr1To19)

      data_merge  = t(rbind(dataExpNorm,dataCnv,dataSom,dataSnp))
    }

  return(list(data_merge=data_merge,
              genename = genename))
}

subStr1To19 = function(x){ substr(x,1,19)}

normalize = function(x){
  x = sapply(x,as.numeric)
  return((x - mean(x)) / sd(x))
}

toRange01 = function(x){
  x = sapply(x,as.numeric)
  y = (x - min(x) + 0.01 * min(x))/(max(x) + 0.01 * max(x) - min(x))
  z = sapply(y,function(xx) min(xx,1))
  return(z)
}

sigmoid  = function(x){
  #function to normalize and transform gene expression using sigmoid function
  x = sapply(x,as.numeric)
  return(exp(x)/(1+exp(x)))
}

logit = function(x){
  x = sapply(x,as.numeric)
  return(log(x /(1 - x)))
}

getSVD    = function(kcval){
  kcval_svdObj  = svd(kcval)
  temp          = cumsum(kcval_svdObj$d)
  cutoff_k      = min(which(temp > 0.8 * sum(kcval_svdObj$d),arr.ind=T))
  return(kcval_svdObj$u[,1:cutoff_k])
}

bestK = function(kcval_svd){
  Ks          = 1: floor(nrow(kcval_svd)/3)
  aic         = (nrow(kcval_svd) - 1) * sum(apply(kcval_svd,2,var))
  aicdiff     = (nrow(kcval_svd) - 2) * sum(apply(kcval_svd,2,var))
  lambda = 0.95
  for (i in Ks) {
    set.seed(12345)
    km      = kmeans(kcval_svd, centers=i)
    #AIC lambda = M (dimension)
    aic[i]  = km$tot.withinss + lambda * i
  }
  
  k_best  = which.min(aic)
  
  return(k_best=k_best)
}

##----------------------------
#
getGroup = function(kcval_svd,kcval){
  k_best                = bestK(kcval_svd)
  colnames(kcval_svd)   = colnames(kcval)[1:ncol(kcval_svd)]
  row.names(kcval_svd)  = colnames(kcval)
  group                 = kmeans(x=kcval_svd,centers=k_best)$cluster
  return(unlist(group))
}

##----------------------------
#
fitCnv = function(X,y,plotflag){
  linFit = lm(y ~ 1 + X)
  beta     = coef(linFit,s=min(ridgefit$lambda))
  # residual = y - beta[1]- X * beta[2]
  residual = linFit$residuals 
  RSS      = sum((residual)^2)       

  if(plotflag == 1){
    plot(y=y,x=X,xlab="cnv",ylab="normalized expression")
    curve (beta[1] + beta[2]*x, add=TRUE,col="blue")
    curve (cbind(1,x) %*% beta, add=TRUE,col="red")
  }
  
  return(list(beta   =beta[-1],
              RSS    =RSS,
              residuals = residual))
}

##cohen's kappe- svd - k-means grouping of variable
groupVars = function(dataSnp,dataCnv,plotflag){
  library(irr)
  kcval = matrix(nrow=nrow(dataSnp), ncol= nrow(dataSnp))
  ###--TODO make it a parallel computing
  for ( i in 1:nrow(dataSnp))  
    for ( j in i:nrow(dataSnp))    
    {
      kp2 = kappa2(ratings= t(rbind(dataSnp[i,],dataSnp[j,])) )
      kcval[i,j]  = kcval[j,i] = kp2$value
    }
  
  colnames(kcval) = row.names(dataSnp)
  rownames(kcval) = row.names(dataSnp)
  kcval_svd       = getSVD(kcval)
  group           = getGroup(kcval_svd,kcval)
  if (nrow(dataCnv) ==1){
    group           = c(cnv = 1,group + 1)
  }
  if (plotflag == 1){
    # pdf(paste(outputcoeff, genename,".pdf",sep=""))
    require(gplots)
    # plot(density(dataExp),main = "Gene Expression Density")
    # plot(density(dataExpSg), main = "Gene Expression Density after normalization")
    heatmap.2(kcval,trace="none",
              main= "kappa similarity matrix")
    heatmap.2(kcval_svd, trace="none",
              main="kappa svd top k dim")
    
  }
  return(group)
}


regfit      = function(X, y, group, plotflag=0){
  #Desp:    function to perform regression 
  #input:   datamatrix with first column as response, other as variables
  #output:  coefficients beta, residual, r_square
  #usage:   function internal called by regfitPermu
  if (setequal(colnames(X), names(group))){
    #search for alpha:
    rss   = NULL
    a     = seq(0.01,1,0.1)
    for (alpha in a){
      grpfit            = grpreg(X,y,group=group,penalty="grLasso",
                                 alpha = alpha,eps=0.005,max.iter=10000)
      grpfit_selected   = select(grpfit,criterion="AIC")
      rss               = c(rss,sum(y - grpfit_selected$beta[1] -
                            X %*% as.matrix(grpfit_selected$beta[-1])))
    }
    alpha = a[which.min(rss)[1]]
    #fit model
    grpfit              = grpreg(X, y, group=group,penalty="grLasso",
                                 alpha = alpha,eps=0.005,max.iter=5000)
    grpfit_selected     = select(grpfit,criterion="AIC")
    
    beta                = grpfit_selected$beta
    residual            = y - beta[1] - X %*% as.matrix(beta[-1])
    residual_sum_squre  = sum((residual)^2)
    reg_sum_square      = sum((X %*% as.matrix(beta[-1]) + beta[1] - mean(y))^2)
    total_sum_square    = reg_sum_square + residual_sum_squre
    r_square            = reg_sum_square / total_sum_square  
    
    if (plotflag == 1){
      plot(grpfit)
      plot(grpfit_selected$IC,xlab="Index of Lambda",ylab="IC",
           sub=paste("Min Lambda", grpfit_selected$lambda))
      qqplot(x=residual,
             y=X %*% as.matrix(beta[-1]) + beta[1],ylab="predicted",
             main="fitted model residual distribution")
      plot(density(residual),
           sub=paste("alpha=",alpha,"lambda=",grpfit_selected$lambda),
           main="Fitted Residual Density")
    }
    return(list(beta    = beta,
                residual          = residual, 
                r_square          = r_square,
                RSS               = residual_sum_squre))
  }
  else { print("ERROR! group not match data!")
  }
}

regfitPermu = function(data_merge,group,nperm, plotflag=0) {
  #Desp: function to calculate the significance of one fitting
  #input: data_merge, group, and number of permuatation, plotflag 1/0
  #output: p-value, plot of p-value significance depends on plotflag
  y               =  sigmoid(data_merge[,1] )
  X               =  as.matrix(data_merge[,-1])
  fit_grpreg      =  regfit(X, y, group, plotflag)
  
  #permutation
  count_perm      = 0 
  r_square_perm   = 0.0
  while (count_perm < as.numeric(nperm) ) {
    y_permu       =  sigmoid(normalize(y + sample(fit_grpreg$residual,size=length(y),replace=TRUE)))
    #y_permu       =  sigmoid(normalize(logit(y) + sample(fit_grpreg$residual,size=length(y),replace=TRUE)))
    if (plotflag == 1){
      if(count_perm == 0){
        plot(density(y_permu),
             ylim = c(0, 1.2 * max(density(y)[[2]])),
             main="permutated Y density",col="gray")
        #lines(density(y_permu),col = "gray")
      }
      else{lines(density(y_permu),col="gray")}
    }
    # residual_perm             = y_permu - X %*% as.matrix(fit_grpreg$beta)
    # residual_sum_squre_perm   = sum((residual_perm)^2)
    # reg_sum_square_perm       = sum((X %*% as.matrix(fit_grpreg$beta) - mean(y_permu))^2)
    # total_sum_square_perm     = reg_sum_square_perm + residual_sum_squre_perm
    # temp                      = 1 - reg_sum_square_perm / total_sum_square_perm 
    # r_square_perm             =   c(r_square_perm, temp) 
    r_square_perm =   c(r_square_perm, regfit(X, y_permu, group)$r_square) 
    count_perm                =   count_perm  + 1
    if(count_perm > 100 && (count_perm %% 100) == 0) {
      print(paste("Permutation", count_perm))}
  }
  if(plotflag == 1) lines(density(y),col="red")
  
  r_square_perm   = r_square_perm[-1]
  pvalue          = max(1/nperm, 
                        length(r_square_perm[which(r_square_perm > fit_grpreg$r_square)]) / nperm)
  
  if (plotflag == 1){
    plot(density(r_square_perm),
         xlab = "" ,
         sub  = paste("pvalue =",pvalue, "\n nperm =",nperm),
         main = "Permutation results")
    abline(v=fit_grpreg$r_square,col="red")
  }
  return(list(pvalue   = pvalue,
              RSS      = fit_grpreg$RSS,
              fit_residual = fit_grpreg$residual,
              fit_r2   = fit_grpreg$r_square,
              beta     = fit_grpreg$beta,
              permu_r2 = r_square_perm,
              grpFit   = fit_grpreg ))
  
}

myFtest = function(residual1,residual2,p1,p2,n){
  #var.test(residual1,residual2,alternative="g")
  f = ((sum(residual1^2) - sum(residual2^2))/(p2-p1)) /(sum(residual2^2)/(n-p2+0.1))
  pval = pf(f,p2-p1,n-p2,lower.tail=F)
  return(list(f=f,pval=pval))
}


