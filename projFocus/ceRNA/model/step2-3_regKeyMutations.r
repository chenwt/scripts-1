setRootd = function(){
  sysInfo = Sys.info()
  if(sysInfo['sysname']=="Darwin" ){
    print("working from MacOS")
    rootd = "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  }else if(sysInfo['sysnameß']=="Linux" ){
    print("working from Linux")
    rootd = "/ifs/home/c2b2/ac_lab/jh3283/projFocus/"
  }
  return(rootd)
}
rootd = setRootd()
source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

wd          = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/CHEK1-temp/"
figd = paste(rootd,"/DATA/projFocus/report/topDown_02042014/fig/",sep="")
cdt = paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="-")

require(gplots)

trimSmpName = function (c) {
  return(vapply(c,FUN=function(x){substr(x,6,15)},'a'))
}

##---init
gene = tail(unlist(strsplit(wd,split="/|-",perl=T)),2)[1]
setwd(wd)
gene = 'KIF23'
fsmps = "samples.txt"
fexp = "exp.mat"
fsom = "som.mat"
fcnv = "cnv.mat"
fmeth = "meth.mat"

##loading data
rawExp = read.delim2(fexp)
smpExp = colnames(rawExp)
colnames(rawExp) = smpExp
gExp = rownames(rawExp)

rawSom = read.delim2(fsom)
smpSom = trimSmpName(colnames(rawSom))
colnames(rawSom) = smpSom
gSom = rawSom[,1]
gMut = apply(rawSom[,1:4],1,function(x){paste(x,collapse=":")})


rawCnv = read.delim2(fcnv)
smpCnv = trimSmpName(colnames(rawCnv))
colnames(rawCnv) =  smpCnv
gCnv = as.character(unlist(rawCnv[,1]))

# rawSnp = read.delim2(fsnp)
# smpSnp =  trimSmpName(colnames(rawSnp))
# colnames(rawSnp) = smpSnp
# gSnp = rawSnp[,1:2]


rawSmps = vapply(unlist(read.table("samples.txt")),as.character,'a')
# smps = Reduce(intersect,list(rawSmps,smpCnv,smpSom,smpSnp,smpExp))
smps = Reduce(intersect,list(rawSmps,smpCnv,smpSom,smpExp))

#-------- DEBUG naming issue
candReg = 'KIF23'
dataExp = apply(subset(rawExp,subset=rownames(rawExp)==candReg, select=smps),2,as.numeric  )
names(dataExp) = smps

idx = grep(gene,gSom,perl=T)
gMutRegulator = gMut[idx]
dataSom = rawSom[idx,]
dataSom = apply(subset(rawSom[idx,],select=smps),2,as.numeric  )
rownames(dataSom) = gMutRegulator
dataSom = t(dataSom[apply(dataSom!=0,1,any),,drop=FALSE])
# colnames(dataSom)

dataCnv = subset(rawCnv,select=smps)
dataCnv = t(dataCnv[apply(dataCnv!=0,1,any),,drop=FALSE])
colnames(dataCnv) = gCnv


# idx = (grep(gene,gSnp[,1],perl=T))
# dataSnp = apply(subset(rawSnp[idx,],select=smps),2,as.numeric)
# dataSnp = t(dataSnp[apply(dataSnp!=0,1,any),,drop=FALSE])
# colnames(dataSnp) = gSnp

###-----work on mutation
mutPosS = vapply(colnames(dataSom),
                 FUN=function(x){as.numeric(unlist(strsplit(x,":"))[3])},1)
mutPosE = vapply(colnames(dataSom),
                 FUN=function(x){as.numeric(unlist(strsplit(x,":"))[4])},1)
mutLen = vapply(colnames(dataSom),
                FUN=function(x){as.numeric(unlist(strsplit(x,":"))[4]) -
                                  as.numeric(unlist(strsplit(x,":"))[3]) },1)
expSomCor.spm = data.frame(matrix(0,nrow=ncol(dataSom)))

res.tt = apply(dataSom, 2, function(x){
  # t.test for differential expression
  tempidx = which(x!=0,arr.ind=T)
    if (length(tempidx) > 2){
        tempres = t.test(dataExp[tempidx],dataExp[-tempidx],
           paired = F,alternative="two.sided")
#         tempcor = cor(dataExp[tempidx],dataExp[-tempidx], method = "spearman")
        return(c(tempres$statistic,tempres$p.value))}
    else{
        return(c(NA,NA))
    }
  })
##p adj
res.tt[2,!is.na(res.tt[2,])] = p.adjust(na.omit(res.tt[2,]),method="fdr")
res.ttt = t(res.tt);colnames(res.ttt) = c("t.stat","pval")
idx = which(as.numeric(res.ttt[,2])< 0.01,arr.ind=T)
corr = sign(res.ttt[,2]) * - log(res.ttt[,2])
corrCol = val2col(corr,col=colorRampPalette(c("steelblue3","white","maroon3"))(25))

pdf(paste(figd,"corrrPlot_KIF23_mutation1k_ttest",cdt,".pdf",sep=""))
plot(mutPosS,corr,
     type="l",col="black", cex.lab = 1.2, font.lab= 2,
     xlab = "mutation hotspot position",ylab = "sign(t.static) * -log(pval)",
     main = paste("test of differential expressino \n caused by mutation hotspot :",gene))
points(mutPosS,corr, col=corrCol,pch=16) #col="red"
points(mutPosS[idx],corr[idx],col="red",pch=23,cex=1.5,bg = "red") #col="red"
text(mutPosS[idx],corr[idx] * 0.95,
     col="black",font=2,
     labels = freq1k[idx])

dev.off()



####---------spearmand coorelation
require("scales")
expSomCor.spm = data.frame(matrix(0,nrow=ncol(dataSom)))
fisherr2z = function(r){
  r == as.numeric(r)
  if(r==0|abs(r)==1){
    return(r)
  }else{
    return(0.5 * (log(1+r) - log(1-r)))
  }
}
res.spmcor = apply(dataSom, 2, function(x){
  ##spearman correlation version
  tempidx = which(x!=0,arr.ind=T)
  if (length(tempidx) >= 3){
    x[-tempidx] = rnorm(length(x[-tempidx]),
                        mean=min(x[tempidx])/10000,
                        sd=2*sd(x[tempidx]))
    tempres = cor.test(dataExp,x,
                       method="spearman")
    r.fisherr2z = fisherr2z(tempres$estimate)
    return(c(tempres$estimate, r.fisherr2z))}
  else{
    return(c(NA,NA))
  }
})
boxplot(corr[!is.na(corr)],plot=F)$out
freq1k = apply(dataSom,2,function(x){
  length(x[x!=0])
  })

corr = res.spmcor[2,]
res.ttt = t(res.spmcor);colnames(res.ttt) = c("corr","corr.adj")
# idx = which(as.numeric(res.spmcor[,2])< 0.01&abs(res.spmcor[,])>=0.5,arr.ind=T)
idx = which(abs(res.ttt[,2])>=0.4,arr.ind=T)
corrCol = val2col(corr,col=colorRampPalette(c("steelblue3","white","maroon3"))(25))
length(corr)

pdf(paste(figd,"corrrPlot_KIF23_mutation1k_spearman",cdt,".pdf",sep=""))
plot(mutPosS,corr,
     type="l",col="black", cex.lab = 1.2, font.lab= 2,
     xlab = "mutation hotspot position",ylab = "corrlation",
     main = paste("test of correlation for: ",gene,"\nexpression and mutation hotspots"))
points(mutPosS,corr, col=corrCol,pch=16) #col="red"
points(mutPosS[idx],corr[idx],col="red",pch=23,cex=1.5,bg = "red") #col="red"
text(mutPosS[idx],corr[idx] * 0.95,
     col="black",font=2,
     labels = freq1k[idx])
dev.off()


##extract mut infor.
res.ttt[idx,]
mapfile = paste(rootd,"/DATA/projFocus/result/02022014/som/test/KIF23_mut1k_mapping.txt",sep="")
dmutmap =  read.table(mapfile)
head(dmutmap)
match(mutPosS[idx[1]],dmutmap$V3)

#### clustering of associations 




print("grouping mutations")
  require("lsa")
  dDist         = cosine(dataSom)  
  colnames(dDist)=rownames(dDist)=colnames(dataSom)
  vars = colnames(dataSom)
  group         = getGroup(t(dataSom),dDist)
  names(group)  = vars
  cntsample = nrow(dataSom) ; cntMut = ncol(dataSom) ;
  require(grpreg)
  require(biclust)
  X = as.matrix(dataSom)
  y = as.matrix(dataExp)
  mysom.gl = regfit(X,y,group)
  reg.grplasso = grpreg(X,y,group=group,penalty="grLasso",
       alpha = 1,eps=0.0005,max.iter=10000)
  par(mfrow=c(1,1))
  plot(reg.grplasso)
  summary(reg.grplasso)
  select(reg.grplasso$lambda)
  rownames(reg.grplasso$beta)[which(reg.grplasso$beta[,1] > 0)]
#    distsvd       = getSVD(dDist)
#   group         = getGroupKmeans(dDist)
#   install.packages("kernlab")
  require(kernlab)
  x = tsne(dataSom,k=3)
  reg.ksvm = ksvm(t(dataSom),scale=F)
  par(mfrow=c(1,2))

  dataSom = dataSom
  runKsvm <- function (dataSom) {
    reg.ksvm = ksvm(t(dataSom),kernel="rbfdot", scale=F)
    candMuts = cbind(colnames(dataSom)[reg.ksvm@alphaindex],
                     reg.ksvm@alphaindex,
                     reg.ksvm@alpha)
    myfml = as.formula(dataExp ~ t(t(dataSom[,reg.ksvm@alphaindex]) * reg.ksvm@alpha))
    reg.lm = lm(myfml)
    summary(reg.lm)
    par(mfrow=c(2,2))
    plot(reg.lm)
    predict(reg.ksvm,)
    table(ifelse(as.numeric(candMuts[,3])==1,"svMut","candsvMut"))
    data.tsne = tsne(t(dataSom),k=3)
    mycol = colorRampPalette(c("blue","red"))(256)
    svcolor = val2col(reg.ksvm@alpha,col=mycol)
    allPcolor = rep("gray",times=nrow(data.tsne))
    allPcolor[reg.ksvm@alphaindex] = svcolor
    require("rgl")
    plot3d( x=data.tsne[,1],y=data.tsne[,2],z=data.tsne[,3], 
            type="s",size=0.8,box=F,
           xlab="",ylab="",zlab="",
           col=allPcolor)
    scatter3D(x=data.tsne[,1],y=data.tsne[,2],z=data.tsne[,3],pch=20,        
              col=allPcolor)
    
    
    
    for (i in 1:length(svcolor)){
      idx = reg.ksvm@alphaindex[i];
      points(datadist.tsne[idx,],col=svcolor[i])
    }
  }

  require("misc3d")ßß
  xyz.coords(x=data.tsne[,1],y=data.tsne[,2],z=data.tsne[,3])
  require("plot3D")
par(mfrow=c(1,1))
  
  scatter3D(x=datadist.tsne[,1],y=datadist.tsne[,2],z=datadist.tsne[,3])

  reg.ksvm = ksvm(dDist,kernel="rbfdot",scale=F)
  reg.ksvm@alpha
  reg.ksvm@alphaindex
  gMutRegulator[reg.ksvm@alphaindex]
  datadist.tsne = tsne(dDist,k=3)
  install.packages("rgl", dependencies = TRUE)
 
  lines3d(x=data.tsne[,1],y=data.tsne[,2],z=data.tsne[,3])
  surface3d(x=data.tsne[,1],y=data.tsne[,2],z=data.tsne[,3])

plot3d(x=datadist.tsne[,1],y=datadist.tsne[,2],z=datadist.tsne[,3])
lines3d(x=datadist.tsne[,1],y=datadist.tsne[,2],z=datadist.tsne[,3])
surface3d(x=datadist.tsne[,1],y=datadist.tsne[,2],z=data.tsne[,3])

  plot(datadist.tsne,type="p")
  mycol = colorRampPalette(c("white","blue"))(256)
  svcolor = val2col(reg.ksvm@alpha,col=mycol)
  for (i in 1:length(svcolor)){
    idx = reg.ksvm@alphaindex[i];
    points(datadist.tsne[idx,],col=svcolor[i])
  }
  lines(x,predict(reg.ksvm,x),col="red")

  library(e1071)
  group.db = dbscan(dDist,eps=0.002,method="dist")
  group.db$cluster
  cntsample = nrow(dataExp) ; cntReg = ncol(dataExp) -1;
require(grpreg)
myexp.gl = regfit(as.matrix(dataExp[,-1]),as.matrix(dataExp[,1]),group)

  
  
#   if(plotflag == 1){
#     image(z=dDist,col=mycol,main="image: similarity som")
#     heatmap.2(dDist,trace="none",col=mycol,dendrogram="none", main="som similarity")
#     image(z=distsvd,col=mycol,main="image: similarity som after svd")
#     row.names(distsvd) = rownames(dDist)
#     heatmap.2(distsvd,col=mycol,trace='none',dendrogram="none", main="heatmap: som similarity after svd",labCol=" ")
#   }
#   
  cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp); cntcnv = ncol(dataCnv); cntsom = ncol(dataSom)
  group     = c(som = rep(1,cntsom), cnv = rep(2,cntcnv), group + 2)
  
  data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp+cntcnv+cntsom))
  colnames(data_merge) = c('exp',colnames(dataSnp),colnames(dataSom),colnames(dataCnv))
  row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames
  
  data_merge[,1] = dataExpß
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


