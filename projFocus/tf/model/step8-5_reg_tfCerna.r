##status: under work

rm(list=ls())

setRootd  = function(){
  sysInfo = Sys.info()
  if(sysInfo['sysname']=="Darwin" ){
    print("working from MacOS")
    rootd = "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  }else if(sysInfo['sysname']=="Linux" ){
    print("working from Linux")
    rootd = "/ifs/home/c2b2/ac_lab/jh3283/"
  }
  return(rootd)
}
rootd     = setRootd()

figd = paste(rootd, "/DATA/projFocus/report/Aug2014/fig",sep="")

source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

###----func-------
formatTestOut = function(test_aov){
  Res.DF_cnv_all = paste(test_aov$Res.Df, collapse="/")
  RSS_cnv_all = paste(round(test_aov$RSS,digits=2),collapse="/")  
  result = c(Res.DF_cnv_all= Res.DF_cnv_all, RSS_cnv_all= RSS_cnv_all, round(unlist(data.frame(test_aov)[2,3:6]),5))
  return(result)
}

orderSampleLable = function(dataExp){
  dataExp$sample <- as.character(dataExp$sample)
  #Then turn it back into an ordered factor
  dataExp$sample <- factor(dataExp$sample, levels=unique(dataExp$sample))
  return(dataExp)
}

noseHeatmap = function(inputdata){
  data = inputdata
  allFeatures = colnames(data)
  data$sample = rownames(data)
  ftr.cRegExp = allFeatures[grep("cRegExp", allFeatures)]
  ftr.cRegCNV = allFeatures[grep("cRegCNV", allFeatures)]
  ftr.cRegTFmut = allFeatures[grep("cRegTFmut", allFeatures)]
  data[,c('cTarExp', ftr.cRegExp, ftr.cRegCNV, ftr.cRegTFmut)] = apply(data[,c('cTarExp', ftr.cRegExp, ftr.cRegCNV, ftr.cRegTFmut)],2,as.numeric)
  
  library(reshape); library(ggplot2) ;library(scales);library(RCurl);  library (grid)
  
  ####---underDevelopment
  
  ### ceRNA target expression
  data.cTarExp = data[,c('sample', 'cTarExp')]
  data.cTarExp$variable = rep("cTarExp     ",NROW(data)); colnames(data.cTarExp) = c("sample", 'value', 'variable')
  data.cTarExp$value <- as.numeric(as.character(data.cTarExp$value))
  data.cTarExp <- orderSampleLable(data.cTarExp)
  data.cTarExp$rescale <- rescale(data.cTarExp$value,to=c(-1,1))
  head(data.cTarExp)
  
  ### ceRNA regulator expression
  data.cRegExp = data[,c('sample', ftr.cRegExp)]
  data.cRegExp[, ftr.cRegExp] = apply(data.cRegExp[,ftr.cRegExp],2, function(x){rescale(x,to=c(-1,1))})
  data.cRegExp.m <- melt(data.cRegExp, id.vars='sample')
  head(data.cRegExp.m) 
  data.cRegExp.m <- orderSampleLable(data.cRegExp.m)
  head(data.cRegExp.m)

  
  ## ceRNA regulator's CNV 
  data.cRegCNV = data[,c('sample', ftr.cRegCNV)]
  data.cRegCNV = melt(data.cRegCNV, id.vars='sample')
  data.cRegCNV$value <- factor(sign(as.numeric(as.character(data.cRegCNV$value))), levels = c(-1,0,1))
  data.cRegCNV <- orderSampleLable(data.cRegCNV)
  head(data.cRegCNV)
  
  ## ceRNA regulator's TF mutation profile (collapsed)
  data.cRegTFmut  = data[,c('sample',ftr.cRegTFmut)]; data.cRegTFmut[,ftr.cRegTFmut] <- sign(data.cRegTFmut[,ftr.cRegTFmut])
  data.cRegTFmut.m = melt(data.cRegTFmut,id.vars='sample')
  data.cRegTFmut.m$value <- factor(data.cRegTFmut.m$value, levels=c(0,1))
  data.cRegTFmut.m <- orderSampleLable(data.cRegTFmut.m)
  
  head(data.cRegTFmut.m)
  
  ### ploting
    base_size <- 6
    p <- ggplot(data.cTarExp, aes(variable, sample) )+ geom_tile(aes(fill = rescale), color = 'white') 
    p1 <- p + scale_fill_gradient2(low='green',mid='white',high='orange') + 
          theme_grey(base_size = base_size) + theme(axis.title=element_blank(), 
#               axis.text.y = element_blank(), axis.text.x = element_blank(),            
              legend.position="none",
              axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
              axis.title.x = element_text(face="bold", colour="#990000", size=8), axis.text.x  = element_text(angle=90, vjust=0.5, size=8))


    p <- ggplot(data.cRegExp.m, aes(variable, sample)) +   geom_tile(aes(fill = value), color = 'white')
    p2 <- p + scale_fill_gradient2(low='green',mid='white',high='orange') + 
            theme_grey(base_size = base_size) +  theme(axis.title=element_blank(), 
#             axis.text.y = element_blank(),axis.text.x = element_blank(), 
            legend.position="none" ,
            axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=8), axis.text.x  = element_text(angle=90, vjust=0.5, size=8))
        
    p <- ggplot(data.cRegCNV, aes(variable, sample)) +   geom_tile(aes(fill = value), color = 'white')
    p3 <- p + scale_fill_manual(values=c("lightblue","white",colors()[119])) +  
#           theme_grey(base_size = base_size) +   
            theme(axis.title=element_blank(), 
#                   axis.text.y = element_blank(),
                    legend.position="none", 
                    axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=8), axis.text.x  = element_text(angle=90, vjust=0.5, size=8))
      
    p <- ggplot(data.cRegTFmut.m, aes(variable, sample)) +  geom_tile(aes(fill = value), colour =   "white") 
    p4 <- p + scale_fill_manual(values=c('white',"red")) + 
            theme(axis.title=element_blank(),  
      #             axis.text.y = element_blank(), 
                    legend.position="none",
            axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=8), axis.text.x  = element_text(angle=90, vjust=0.5, size=8))
    
     
  len1 = length(unique(data.cTarExp$variable)); len2 = length(unique(data.cRegExp.m$variable)); len3= length(unique(data.cRegCNV$variable)); len4=length(unique(data.cRegTFmut.m$variable)); 
  lensum = len1 + len2 + len3 + len4
  pdf(paste(figd,"/heatmap_diagnois_regCerna_Aug17_",tgene, ".pdf",sep=""),width=15,height=20)
#   layt<-grid.layout(nrow=1,ncol=4,widths=c(len1/lensum,len2/lensum,len3/lensum,len4/lensum))
  
  layt<-grid.layout(nrow=1,ncol=4,widths=c(1/10,3/10,3/10,3/10))
  
  #   grid.show.layout(layt)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(p1,  vp = viewport(layout.pos.row=1,layout.pos.col=1))
  print(p2, vp = viewport(layout.pos.row=1,layout.pos.col=2))
  print(p3, vp = viewport(layout.pos.row=1,layout.pos.col=3))
  print(p4, vp = viewport(layout.pos.row=1,layout.pos.col=4))
  dev.off()

}

tarTFmutReg = function(tgene, plot = "False"){
  inputdata= paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/",tgene,".tgTFreg.input",sep="")
  inputdata2 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/tfGreedy_",tgene,".temp",sep="")
  inputdata3 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/tfGreedy_act_",tgene,"_1000_1000_mutSampleVector",sep="")
  
  mutSmp = read.table(inputdata3, header=T,stringsAsFactors=F); smpMerg = colnames(mutSmp)[which(mutSmp[1,]>0)]
  mutD = read.table(inputdata2,header=T,stringsAsFactors=F); dataExp = t(mutD[1,]); mutD = data.frame(t(mutD[-1,] )); 
  allsmp1 = rownames(mutD); mutD[setdiff(allsmp1,smpMerg),] <- 0
  length(which(rowSums(mutD)>0))
  
  datacnv = t(read.table(inputdata, header=T,stringsAsFactors=F))
  head(datacnv)
  allsmp2 = rownames(datacnv)
  allsmp = intersect(allsmp1, allsmp2); datacnv <- datacnv[allsmp,]
  
  data = data.frame(cbind(dataExp, mutD)); tgene = colnames(data)[1]; allregs = colnames(data)[-1]; colnames(data) <- c('exp',unlist(as.character(allregs)))
  data$cnv  <- 0; 
  data[allsmp,'cnv'] <-  unlist(as.numeric(datacnv[,2]))
  
  data[setdiff(allsmp,smpMerg),allregs] <-0
  
  ftest = anova(lm(exp~cnv,data),lm(exp~.,data[,c('exp','cnv',allregs)]))
  res = c(tgene = tgene,  formatTestOut(ftest),allregs = paste(allregs, collapse=","))
  if (plot == 'True'){noseHeatmap(data)}
  return(res)
}

testNormality = function(tgene){
  ### exp-cnv--mut(wrong)
  inputdata= paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/",tgene,".tgTFreg.input",sep="")
  ### exp-mut(beforeGreedy)
  inputdata2 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/tfGreedy_",tgene,".temp",sep="")
  ### target-collapsedTFMut(afterGreedy)
  inputdata3 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/tfGreedy_act_",tgene,"_1000_1000_mutSampleVector",sep="")
  print(tgene)
  if (!(file.exists(inputdata) && file.exists(inputdata2) && file.exists(inputdata3))) {
    return(c('NA','NA','NA'))
  }
  mutSmp = read.table(inputdata3, header=T,stringsAsFactors=F); smpMerg = colnames(mutSmp)[which(mutSmp[1,]>0)]
  mutD = read.table(inputdata2,header=T,stringsAsFactors=F); dataExp = t(mutD[1,]); mutD = data.frame(t(mutD[-1,] )); 
  allsmp1 = rownames(mutD); mutD[setdiff(allsmp1,smpMerg),] <- 0
  length(which(rowSums(mutD)>0))
  dataExp = rescale(dataExp,to=c(1,10))
  res1 = c(tgene = tgene, before=shapiro.test(dataExp)$p.value, log2=shapiro.test(log2(dataExp))$p.value)
  #     return(res1)
  datacnv = t(read.table(inputdata, header=T,stringsAsFactors=F))
  head(datacnv)
  allsmp2 = rownames(datacnv)
  allsmp = intersect(allsmp1, allsmp2); datacnv <- datacnv[allsmp,]
  
  data = data.frame(cbind(dataExp, mutD)); tgene = colnames(data)[1]; allregs = colnames(data)[-1]; colnames(data) <- c('exp',unlist(as.character(allregs)))
  data$cnv  <- 0; 
  data[allsmp,'cnv'] <-  unlist(as.numeric(datacnv[,2]))
  
  data[setdiff(allsmp,smpMerg),allregs] <-0
  
  data$exp <- log2(data$exp)
  ftest = anova(lm(exp~cnv,data),lm(exp~.,data[,c('exp','cnv',allregs)]))
  res2 = c(formatTestOut(ftest),allregs = paste(allregs, collapse=","))
  
  return(c(res1,res2))
}


###----main-------
require(scales)

# gfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfCerna/test/input.gene"
# tgenelist = unlist(read.table(gfile,stringsAsFactors=F,header=T)[,1])
# tgene = tgenelist[1]
# 
# cnvfile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv/brca_cnv_l3_tumor_Mar-23-2014.matrix.uniq.matrix"
# expfile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
# mutActfile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix"
# mutOptfileCollapse="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/runJuly27/summary/all_targets_mutSampleVector"
# cernafile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_ceRNA_network.txt"

inputfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfCerna/test/input_PGM2L1"
data = read.csv2(inputfile, sep="\t", header=T,stringsAsFactors=F)
data = data[order(as.numeric(as.character(data$cTarExp))),]
allSmp = data$X; rownames(data) <- allSmp
data = data[,colSums(data!=0)!=0]


allFeatures = colnames(data)
ftr.cTarExp = allFeatures[grep("cTarExp", allFeatures)]
ftr.cRegExp = allFeatures[grep("cRegExp", allFeatures)]
ftr.cRegCNV = allFeatures[grep("cRegCNV", allFeatures)]
ftr.cRegTFmut = allFeatures[grep("cRegTFmut", allFeatures)]

data[,c(ftr.cRegTFmut, ftr.cRegCNV)] <- sign(apply(data[,c(ftr.cRegTFmut, ftr.cRegCNV)],2, as.numeric))
data[,c(ftr.cTarExp, ftr.cRegExp)] <- apply(data[,c(ftr.cTarExp, ftr.cRegExp)], 2, as.numeric)



### plot
# noseHeatmap(data)

require(glmnet)

### elastic net, cRegCNV V.S cRegTFmut + cRegCNV 

data.mod1 = as.matrix(data[,c(ftr.cTarExp, ftr.cRegCNV)])
data.mod2 = as.matrix(data[,c(ftr.cTarExp, ftr.cRegCNV, ftr.cRegTFmut)])

head(data.mod1)
mod1 <- cv.glmnet(x=data.mod1[,-1],y=data.mod1[,1], family='gaussian',alpha=0, nfolds=10)
mod1.lambdamin  = mod1$lambda.min; idx.lambda1 = which(mod1$lambda==mod1.lambdamin)

idx.exlude = which(mod1$glmnet.fit$beta[,idx.lambda1] !=0)

names(which(mod1$glmnet.fit$beta[,idx.lambda1] !=0))

head(data.mod2)
mod2 <- glmnet(x= data.mod2[,-1], y=data.mod2[,1], family='gaussian',alpha=0.95, exclude=idx.exlude)
mod2 <- cv.glmnet(x= data.mod2[,-1], y=data.mod2[,1], family='gaussian',alpha=0, nfolds=10)
plot(mod2)

mod2.lambdamin  = mod2$lambda.min; idx.lambda2 = which(mod2$lambda==mod2.lambdamin)
names(which(mod2$glmnet.fit$beta[,idx.lambda2] !=0))

getRss = function(mod, data.mod){
  beta = mod$glmnet.fit$beta[,which(mod$lambda == mod$lambda.min)]
  y.hat = rowSums(beta * data.mod[,-1])
  y = data.mod[,1]
  rss = sum((y-y.hat)^2)
  return(rss)
}

rss1 = getRss(mod1,data.mod1)
rss2 = getRss(mod2,data.mod2)
rss1
rss2

idx.minlambda = which(mod2$lambda == mod2$lambda.min); p2 = mod2$glmnet.fit$df[idx.minlambda]

mod1 = glm(cTarExp~., data= data.frame(data.mod1),family='gaussian')
mod2 = glm(cTarExp~., data= data.frame(data.mod2),family='gaussian')

require(mixtools)
require(e1071)
cost = 1000; gamma = 0.001
testReg1 = svm(cTarExp ~ ., data = data.frame(data.mod1), type="eps-regression",  kernel = "polynomial", cost = cost, gamma = gamma)
testReg2 = svm(cTarExp ~ ., data = data.frame(data.mod2), type="eps-regression",  kernel = "polynomial", cost = cost, gamma = gamma)

rss1 = sum(testReg1$residuals) ^2
rss2 = sum(testReg2$residuals)^2

p1 = NCOL(data.mod1)
p2 = NCOL(data.mod2)

n = NROW(data)

f = ((rss1-rss2)/(p2-p1))/(rss2/(n-p2))

##alternative hypo: mod2 is bettern than mod1, rss1 is smaller than rss2
f.pval = ifelse(f>0,pf(f, p2-p1, n-p2, lower.tail=T), 1)

}

