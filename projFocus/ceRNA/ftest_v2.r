#!/usr/bin/Rscript
#J.HE
#input: gene_exp; gene_cnv; gene_som; gene_snp;
#output: file of gene snp contribution to expression significance
#Descripiton: This is to use ridge linear fitting mode to fitting cnv+somatic mutation, and cnv+somatic + snp,
# 			  calculating the significance of SNPs' contribution for gene expression
# Path: /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test
#TODO: add INDEL data handler


sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/fTest/test")
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("loaded from Linux")
  setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/fTest/test")
}

args = getArgs()
usage = "-Usage: Rscript ftest.r --exp exp.mat --snp snp.mat --cnv cnv.mat --som somaticMutation.mat --type 1/2/3[1: snp,cnv,som;2:snp,cnv; 3:snp,som;]
Example: "
if(length(args) < 5 || is.null(args)){
	print(usage)
	stop("Input parameter error!")
}

setwd(system("pwd",intern=T))
cwd         = getwd()
inputsnp    = args['snp']
inputexp    = args['exp']
inputcnv    = args['cnv']
inputsom    = args['som']
inputtype   = args['type'] #1,2,3
outputfile  = args['out']
output      = paste(cwd,"/",outputfile, sep="")
print(paste("current working directory",getwd()))

###-----usefule functions----
formatColname = function(d){
 vapply(colnames(d),FUN=subStr1To19,'a')
}

formatData = function(inputFile,t='exp'){
  switch(t,'exp'={
    data = read.table(inputFile,header=T)
    geneName = as.character(data[1,1])
    row.names(data) = geneName
    colnames(data) = formatColname(data)
    data = data[,-c(1:4)]
  },'snp'={
    data = read.table(inputFile,header=T)
    row.names(data) = as.character(data$snpname)
    colnames(data) = formatColname(data)
    data = data[,-c(1,2)]   
  },'cnv'={
    data = read.table(inputFile,header=T)[,-c(1,2)]
    row.names(data) = rep('cnv',nrow(data))
    colnames(data) = formatColname(data)
  },'som'={
    data = read.table(inputFile,header=T)[,-c(1,2)]
    row.names(data) = rep('som', nrow(data))
    colnames(data) = formatColname(data)
  },'indel'={
    ##under development
  })
  return(t(data))
} 

###-----functions---end-----



###------------test-start-------
# cwd         = getwd()
# inputsnp    =  "ESR1_brca_GWASCataLogGene_snp_KWtest.mat.anno.adjPass_1.0.mat"
# inputexp    =  "ESR1_brca_exp_l3_731_DEG.mat.singleTSS.anno"
# inputcnv    =  "ESR1_brca_gene_DEG_cnv_731.mat"
# inputsom    =  "ESR1_brca_somForDeg.mat"
# inputtype    = 1
# output = paste(cwd,"/fTest_pval.txt",sep="")
# print(paste("current working directory",getwd()))
###------------test-end--------

#load data
dataExp= formatData(inputexp,t='exp')
geneName=colnames(dataExp)
dataSnp = formatData(inputsnp,t='snp')

require(MASS)
switch(inputtype,'1'={
  print("running exp=snp+cnv+som model...")
  #-----model-exp=snp+cnv+som    
  dataSom = formatData(inputsom,t='som')
  dataCnv = formatData(inputcnv,t='cnv')
  
  cntexp= ncol(dataExp); cntsnp = ncol(dataSnp); cntcnv=ncol(dataCnv); cntsom = ncol(dataSom)
  cntsample = nrow(dataExp)
  
  data = as.data.frame(matrix(NA,ncol=(cntexp +cntsnp  + cntcnv + cntsom ), nrow=cntsample))
  row.names(data) = row.names(dataExp)
  colnames(data) = c("exp",colnames(dataSnp), colnames(dataSom), colnames(dataCnv))
  ifelse(cntexp > 1, (data[,1:cntexp]=dataExp), (data[,1]=dataExp))
  data[,(cntexp+ 1):(cntexp + cntsnp) ]= dataSnp
  data[,(cntexp+cntsnp+1):(cntexp+cntsnp+cntsom)] = dataSom
  data[,(cntexp+cntsnp+cntsom + 1) :(cntexp+cntsnp+cntsom+cntcnv) ] = dataCnv
  
  fit1 = lm(exp ~ 1 + som +cnv, data=data)
  fit2 = lm(exp ~ 1 + ., data=data)
  fresult = anova(fit1,fit2)
  names(fresult) = c("ResDf","RSS","DF","SumOfSquare","F","PrGreaterThanF")
  result = fresult$PrGreaterThanF[2]
},
  '2'= {
    print("running exp=snp+cnv model...")    
    #---model-exp=snp+cnv
    dataCnv = formatData(inputcnv,t='cnv')
    
    cntexp= ncol(dataExp); cntsnp = ncol(dataSnp); cntcnv=ncol(dataCnv)
    cntsample = nrow(dataExp)
    
    data = as.data.frame(matrix(NA,ncol=(cntexp + cntsnp  + cntcnv), nrow=cntsample))
    row.names(data) = row.names(dataExp)
    colnames(data) = c("exp",colnames(dataSnp), colnames(dataCnv))
    ifelse(cntexp > 1, (data[,1:cntexp]=dataExp), (data[,1]=dataExp))
    data[,(cntexp+ 1):(cntexp + cntsnp) ]= dataSnp
    data[,(cntexp+cntsnp+ 1) :(cntexp+cntsnp+cntcnv) ] = dataCnv
    
    fit1 = lm(exp ~ 1 + cnv, data=data)
    fit2 = lm(exp ~ 1 + ., data=data)
    fresult = anova(fit1,fit2)
    names(fresult) = c("ResDf","RSS","DF","SumOfSquare","F","PrGreaterThanF")
    result = fresult$PrGreaterThanF[2]
    
  },
  '3'={
    print("running exp=snp+som model...")
    #---model-exp=snp+som
    dataSom = formatData(inputsom,t='som')
    
    cntexp= ncol(dataExp); cntsnp = ncol(dataSnp); cntsom = ncol(dataSom)
    cntsample = nrow(dataExp)
    
    data = as.data.frame(matrix(NA,ncol=(cntexp +cntsnp  +  cntsom ), nrow=cntsample))
    row.names(data) = row.names(dataExp)
    colnames(data) = c("exp",colnames(dataSnp), colnames(dataSom))
    ifelse(cntexp > 1, (data[,1:cntexp]=dataExp), (data[,1]=dataExp))
    data[,(cntexp+ 1):(cntexp + cntsnp) ]= dataSnp
    data[,(cntexp+cntsnp+1):(cntexp+cntsnp+cntsom)] = dataSom
    
    fit1 = lm(exp ~ 1 + som , data=data)
    fit2 = lm(exp ~ 1 + ., data=data)
    fresult = anova(fit1,fit2)
    names(fresult) = c("ResDf","RSS","DF","SumOfSquare","F","PrGreaterThanF")
    result = fresult$PrGreaterThanF[2] 
  }
  #--add indel in the future
)

names(result) = paste(geneName,inputtype,sep="_")
write.table(as.matrix(result),file=output,col.names=F,quote=F,sep="\t",append=T)
print("#--DONE---")