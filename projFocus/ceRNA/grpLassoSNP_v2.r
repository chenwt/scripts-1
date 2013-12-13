##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file: expression and snp gt for one gene > 
#output: <file: gene: snp with weight 
#Usage: Rscript grplassoSNP.r input
#Description: require package irr, grpreg,
#             Major change: without sigmoid transformation
#TODO:add quality checking script for each file and the content

#getting command line parameters
plotflag    = 0 # "0  for no plotting "
args        = commandArgs(TRUE)
if (is.null(args) || length(args) < 3){
  print("Please provide parameters in right order!")
  print("Usage: Rscript grplassoSNP_v1.r <inputsnp> <inputexp> <inputcnv> <1/0>")
  stop("Stop since no parameters")
}else{
  print(args)
}

setwd(system("pwd",intern=T))
cwd         = getwd()
inputsnp    = args[1]
inputexp    = args[2]
inputcnv    = args[3]
plotflag    = args[4]
outputcoeff = paste(cwd,"/grplasso_coeff_",sep="")
print(paste("current working directory",getwd()))
# print(paste("Current working directory:",system("pwd",intern=T)))

##-----------------------test-##----------------------------
#setwd("/Volumes/ac_lab/jh3283/scripts/projFocus/ceRNA/test/")
# setwd("/Volumes//ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr22/")
#setwd("/Users/jh3283/projFocus")
# inputsnp    = "temp/input_snp_ADM2"
# inputexp    = "temp/input_exp_ADM2"
# inputcnv    = "temp/input_cnv_ADM2"
# # inputsnp    = "input_test_reg_snp.mat"
# # inputexp    = "input_test_reg_exp.mat"
# # inputcnv    = "input_test_reg_cnv.mat"
#outputcoeff = "grplasso_coeff_"
##----------------------------##----------------------------

##-------useful functions-------------------
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

getGroup = function(kcval_svd){
  k_best                = bestK(kcval_svd)
  colnames(kcval_svd)   = colnames(kcval)[1:ncol(kcval_svd)]
  row.names(kcval_svd)  = colnames(kcval)
  group                 = kmeans(x=kcval_svd,centers=k_best)$cluster
  return(unlist(group))
}

##----------------------------##----------------------------
print("loading data")

##load exp data
dataExp             = read.table(inputexp,header =T)
rownames(dataExp)   = dataExp[,1]
genename            = dataExp[,1] 
dataExp             = sapply(dataExp[,-c(1:4)],
                            function(x){as.numeric(as.character(x))})

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
    dataExpSg   = normalize(dataExp)
    data_merge  = t(rbind(dataExpSg,dataCnv,dataSnp))

    ##--------------------------------Grouping_of_variables
    library(irr)
    print("grouping of variables...")
    kcval = matrix(nrow=nrow(dataSnp), ncol= nrow(dataSnp))
    ###--TODO make it a parallel computing
    for ( i in 1:nrow(dataSnp))  
      for ( j in i:nrow(dataSnp))    
      {
        kp2 = kappa2(ratings= t(rbind(dataSnp[i,],dataSnp[j,])) )
        kcval[i,j]  = kcval[j,i] = kp2$value
        #kcp[i,j]   = kcp[j,i]   = kp2$p.value 
      }
    colnames(kcval) = row.names(dataSnp)
    rownames(kcval) = row.names(dataSnp)
    kcval_svd       = getSVD(kcval)
    group           = getGroup(kcval_svd)

    #check whether cnv data present
    if (nrow(dataCnv) ==1){
      group           = c(cnv = 1,group + 1)
    }

    if (plotflag == 1){
      library("gplots",lib.loc="/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/library")
      pdf(paste(outputcoeff, genename,".pdf",sep=""))
      plot(density(dataExp),main = "Gene Expression Density")
      plot(density(dataExpSg), main = "Gene Expression Density after normalization")
      heatmap.2(kcval,trace="none",
                main= "kappa similarity matrix")
      heatmap.2(kcval_svd, trace="none",
                main="kappa similarity matrix after svd top k dim")
    }

    ##########-----------model------------
    ###group lasso 
    #adding cnv variable
    #TODO: calculate cnv only residual, exclude confunding factor from CNV

    ##---------cnv+snp
    require(grpreg)
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
          rss               = c(rss,sum(y - X %*% as.matrix(grpfit_selected$beta[-1])))
        }
        #fit for best alpha
        alpha = a[which.min(rss)[1]]
        grpfit              = grpreg(X, y, group=group,penalty="grLasso",
                                     alpha = alpha,eps=0.005,max.iter=5000)
        grpfit_selected     = select(grpfit,criterion="AIC")
        
        beta                = grpfit_selected$beta[-1]
        residual            = y - X %*% as.matrix(beta)
        residual_sum_squre  = sum((residual)^2)
        reg_sum_square      = sum((X %*% as.matrix(beta) - mean(y))^2)
        total_sum_square    = reg_sum_square + residual_sum_squre
        r_square            = reg_sum_square / total_sum_square  # new one
        
        #plot 
        if (plotflag == 1){
          plot(grpfit)
          plot(grpfit_selected$IC,xlab="Index of Lambda",ylab="IC",
               sub=paste("Min Lambda", grpfit_selected$lambda))
          qqplot(x=residual,
                 y=X %*% as.matrix(beta),ylab="predicted",main="fitted model residual distribution")
          plot(density(residual),
               sub=paste("alpha=",alpha,"lambda=",grpfit_selected$lambda),
               main="Fitted Residual Density")
        }
        #return values
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
      print("Doing permutation......")
      while (count_perm < as.numeric(nperm) ) {
        y_permu       =  sigmoid(normalize(y + sample(fit_grpreg[[2]],size=length(y),replace=TRUE)))
        r_square_perm =   c(r_square_perm, regfit(X, y_permu, group)$r_square) 
        count_perm                =   count_perm  + 1
        if (plotflag == 1){
          if(count_perm == 1){
            plot(density(y_permu),
                 main="permutated Y density",col="gray",type="n")
            lines(density(y_permu),col="gray")
          }
          else{lines(density(y_permu),col="gray")}
        }
        if(count_perm > 100 && (count_perm %% 100) == 0) {
          print(paste("Permutation", count_perm))}
          
      }
      #add actual data density line
      if (plotflag ==1){lines(density(y),col="red")}
      
      r_square_perm   = r_square_perm[-1]
      pvalue          = max(1/nperm, 
                            length(r_square_perm[which(r_square_perm > fit_grpreg$r_square)]) / nperm)

      if (plotflag == 1){
          plot(density(r_square_perm),
               xlab = "",xlim=c(0,1),
               sub  = paste("pvalue =",pvalue, "\n nperm =",nperm),
               main = "Permutation results")
          abline(v=fit_grpreg$r_square,col="red")
      }
        return(list(pvalue = pvalue,
                    RSS    = fit_grpreg$RSS,
                    fit_r2 = fit_grpreg$r_square,
                    beta   = fit_grpreg$beta,
                  permu_r2 = r_square_perm ))
     
    }

    print("Doing regression...")
    nperm     = 1000
    fitpermut = regfitPermu(data_merge,group,nperm, plotflag)
    # system.time(regfitPermu(data_merge,group,nperm=1000))
    # 80s for 1000 permutation

    ##output
    print("writing output...")
    outputcoeff = paste(outputcoeff,genename,".txt",sep="")
    write.table(t(as.matrix(c(paste("gene","RSS","npermu","pvalue",sep=":"),paste(genename,fitpermut$RSS, nperm, fitpermut$pvalue,sep=":")))),
                outputcoeff,
                quote=F,col.names=F,sep="\t",row.names = F)
    write.table(as.matrix(sort(fitpermut$beta)),
                outputcoeff,
                append=T,
                col.names=F,quote=F,sep="\t")
}else {
  write.table(t(as.matrix(c(paste("gene","RSS","npermu","pvalue",sep=":"),paste(genename,'NA', 'NA', 'NA',sep=":")))),
                  outputcoeff,
                  quote=F,col.names=F,sep="\t",row.names = F)
}