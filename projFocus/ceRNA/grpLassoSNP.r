##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file: expression and snp gt for one gene > 
#output: <file: gene: snp with weight 
#Usage: Rscript grplassoSNP.r input
#Description: require package irr, grpreg,gplots
#             Major change:  add funtion to get command line arguments using --arg value pair
#                           formalized model
#                           
#TODO: develop plot option!

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test")
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test")
}

###-----------------functions-----------------
kappaDist = function(dataSnp){
  library(irr)
  cntsnp = ncol(dataSnp)
  cntsample = nrow(dataSnp)
  kcval = matrix(NA, nrow=cntsnp, ncol= cntsnp)
  ###--TODO make it a parallel computing
  for ( i in 1:cntsnp){
    temp = as.data.frame(matrix(1,cntsample, 2))
    temp[,1] = dataSnp[,i]
    for ( j in i:cntsnp)    
    {temp[,2] = dataSnp[,j]
     kp2 = kappa2(ratings= temp)
     kcval[i,j]  = kcval[j,i] = kp2$value
    }
  }
  colnames(kcval) = colnames(dataSnp)
  rownames(kcval) = colnames(dataSnp)
  return(kcval)
}

imageDist = function(dst){
  dim <- ncol(dst)
  image(1:dim, 1:dim, dst, axes = FALSE) 
  axis(1, 1:dim, nba[1:20,1], cex.axis = 0.5)
  axis(2, 1:dim, nba[1:20,1], cex.axis = 0.5)
  for (i in 1:dim)
  {
    for (j in 1:dim)
    {
      txt <- sprintf("%0.1f", dst[i,j])
      text(i, j, txt, cex=0.5)
    }
  }
}
###------------------functions---end---------

#getting command line parameters
args = getArgs()
plotflag = 0
nperm    = 1000
usage = "Usage: Rscript grpLassoSNP.r --exp exp.mat --snp snp.mat --cnv cnv.mat --som som.mat --out outputfile --type 1[2/3,model type] [optional: --plot 1/0[plot/no plot] --nperm 1000]"
example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r "
if(length(args) < 7 || is.null(args)){
  print(paste(usage,example,sep="\n"))
  print(args)
  stop("Input parameter error!")
}

setwd(system("pwd",intern=T))
cwd         = getwd()
inputsnp    = args['snp']
inputexp    = args['exp']
inputcnv    = args['cnv']
inputsom    = args['som']
type        = args['type']
plotflag    = as.integer(args['plot']) #optional 
# nperm       = as.integer(args['nperm'])
output      = paste(cwd,"/",args['out'], sep="")
print(paste("current working directory:",cwd))
# print(paste(inputsnp,inputexp,inputcnv,inputsom,type,plotflag,nperm,output,sep="\n"))
##-----------------------test-##----------------------------
# setwd(system("pwd",intern=T))
# cwd         = getwd()
# inputsnp    = "input_test_reg_snp.mat"
# inputexp    = "ESR1_brca_exp_l3_731_DEG.mat.singleTSS.anno"
# inputcnv    = "ESR1_brca_gene_DEG_cnv_731.mat"
# inputsom    = "ESR1_brca_somForDeg.mat"
# type        = 1
# plotflag    = 1
# nperm       = 100
# output      ="grplasso_coeff_"
# output      = paste(cwd,"/",output, sep="")
# print(paste("current working directory:",cwd))
##----------------------------##----------------------------
##load exp data
dataExp = formatData(inputFile=inputexp,t='exp')
genename = colnames(dataExp)
dataExp = normalize(dataExp)
outtxt = paste(output,genename,".txt",sep="")
outplot = paste(output,genename,".pdf",sep="")
print(outtxt)
print(outplot)
if (plotflag == 1){
     print("plotting...")
     pdf(outplot)
     require(gplots)
     mycol = bluered(256)
     plot(density(dataExp),main = "Gene Expression Density")
    }
dataSnp = formatData(inputFile=inputsnp,t='snp')
cntsnp = ncol(dataSnp)
if (cntsnp > 1){
  print("grouping snps...")
  kcDist          = kappaDist(dataSnp) 
  kcDist_svd      = getSVD(kcDist)
  group           = getGroup(kcDist_svd,kcDist)
  if(plotflag == 1){
    image(z=kcDist,col=mycol,main="image: KC similarity SNP")
    heatmap.2(kcDist,trace="none",col=mycol,dendrogram="none", main="heatmap kc similarity")
    image(z=kcDist_svd,col=mycol,main="image: KC similarity SNP after svd")
    row.names(kcDist_svd) = rownames(kcDist)
    heatmap.2(kcDist_svd,col=mycol,trace='none',dendrogram="none", main="heatmap: KC similarity SNP after svd",labCol=" ")
  }
}else{
  print("only one snp...")
  group   = as.integer(1)
  names(group) = colnames(dataSnp)
}  
switch(type,'1'={
    ##model 1 exp = snp + som + cnv
    print("model exp~snp+som+cnv")
    dataCnv = formatData(inputcnv,t='cnv')
    dataSom = formatData(inputsom,t='som')
    cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp); cntcnv = ncol(dataCnv); cntsom = ncol(dataSom)
    group     = c(som = rep(1,cntsom), cnv = rep(2,cntcnv), group + 2)
    
    data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp+cntcnv+cntsom))
    colnames(data_merge) = c('exp',colnames(dataSnp),colnames(dataSom),colnames(dataCnv))
    row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames
    
    data_merge[,1] = dataExp
    data_merge[,2:(1+cntsnp)] = dataSnp
    data_merge[,(1+cntsnp+1):(1+cntsnp+cntsom)] = dataSom
    data_merge[,(1+cntsnp+cntsom+1):(1+cntsnp+cntsom+cntcnv)] = dataCnv

    },'2'={    
      ##model 2 exp = snp + soｍ 
      print("model exp~snp + som")
      dataSom = formatData(inputsom,t='som')
      cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp); cntsom = ncol(dataSom)
      group     = c(som = rep(1,cntsom), group + 1)
      
      data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp+cntsom))
      colnames(data_merge) = c('exp',colnames(dataSnp),colnames(dataSom))
      row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames
      
      data_merge[,1] = dataExp
      data_merge[,2:(1+cntsnp)] = dataSnp
      data_merge[,(1+cntsnp+1):(1+cntsnp+cntsom)] = dataSom
        
    },'3'={
     ## model 3 exp = snp + cnv 
      print("model exp~snp + cnv")
      dataCnv = formatData(inputcnv,t='cnv')
      cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp); cntcnv = ncol(dataCnv)
      group     = c( cnv = rep(1,cntcnv), group + 1)
      
      data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp+cntcnv))
      colnames(data_merge) = c('exp',colnames(dataSnp),colnames(dataCnv))
      row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames

      data_merge[,1] = dataExp
      data_merge[,2:(1+cntsnp)] = dataSnp
      data_merge[,(1+cntsnp+1):(1+cntsnp+cntcnv)] = dataCnv      
    }, '4'={
      print("model exp ~ snp")
      cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp);
      data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp))
      colnames(data_merge) = c('exp',colnames(dataSnp))
      row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames    
      data_merge[,1] = dataExp
      data_merge[,2:(1+cntsnp)] = dataSnp
    })
  
require(grpreg)
print("Doing regression...")
nperm     = nperm
fitpermut = regfitPermu(data_merge, group, nperm, plotflag=plotflag) 

##output
print("writing output...")
write.table(t(as.matrix(c(paste("gene","RSS","npermu","pvalue",sep=":"),paste(genename,fitpermut$RSS, nperm, fitpermut$pvalue,sep=":")))),
            outtxt,
            quote=F,col.names=F,sep="\t",row.names = F)
write.table(as.matrix(sort(fitpermut$beta[-1])),
            outtxt,
            append=T,
            col.names=F, quote=F,sep="\t")
print(paste("#-----Done",genename))