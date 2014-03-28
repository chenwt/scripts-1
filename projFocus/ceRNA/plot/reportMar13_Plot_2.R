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
fsnp = "snp.mat"
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

