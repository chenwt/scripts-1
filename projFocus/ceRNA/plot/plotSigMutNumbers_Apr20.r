###-------sig target genesf
dfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/run-Apr-1-2014/summary/keyReg_Apr-20-2014.summary"
figD  = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Apr2014/fig/"
dataF = read.table(dfile, sep = "\t", header = T)
dataF = dataF[which(dataF$regSig>0),]
# dataF[order(dataF$r2,decreasing=T),]
r2T  = table(cut(dataF$r2,breaks=seq(0,1,by=0.1)))
pdf(paste(figD,"sigGintPlot_Apr20.pdf",sep=""))
par(mar=c(6,6,2,1),mgp=c(4.5,0.6,0.5))
barplot(cumsum(r2T)/sum(r2T) * 100,space=0.2,width=1.2,col="lightblue", border="gray", font=2,
        ylim = c(0,101),las=2,
        xlab = "r2 value",ylab = "% regulators", main = "R2 value of Gint's Regression")
text(x=2,y=90, label=paste("total sig Gint: ",length(dataF$r2),sep=""),font=2)

par(mar=c(3.5,3.5,2,1),mgp=c(2,0.6,0.5))
hist(dataF$r2,col="lightblue", border="gray", font=2,freq=F,las=2,
        xlab = "r2 value",ylab = "density", main = "R2 value of Gint's Regression")
lines(density(dataF$r2),col="orange", lwd = 2.5)

dev.off()


###------Total mutated Regulator
dfile = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/keyReg_Apr-20-2014.summary.driverRegs.list.mut.matrix.mutSampleCount"
figD  = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Apr2014/fig/"

data = read.table(dfile, sep="\t",header=F)
dT   = table(data$V2)
dpercentT = round(dT/sum(dT) * 100,3)
pdf(paste(figD,"/barplot_totalMutatedReg_Apr21-2014.pdf",sep=""))
barplot(dpercentT, 
     xlab = "Recurrence", ylab = "%", main = "Total Mutated Regulators", font = 2,
     col = "lightblue", border="gray")
text(seq_len(length(dT)) * 1.2 -  0.5 ,vapply(dpercentT, FUN=function(x){max(5,x*0.5)}, 1),
     col="black",font=2,cex=1.5,
     labels = dT)
dev.off()

###-------Mutated Regulators for each cancer gene
require("gridExtra")
dfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/mutRegSamp.numbers.Apr172014_sig"
figD  = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Apr2014/fig/"

data  = read.table(dfile, sep="\t", header = T)
# pdf(paste(figD,"mutatedSigRegulator_Apri20-2014.pdf",sep=""),width=10)

# hist(data$gintSmpCnt, breaks=20,freq=F,
#      xlab = "Number of Gint samples", ylab = "Frequency", main = "Mutated Gint Samples",
#      col = "lightblue", border="gray")
# lines(density(as.numeric(data$gintSmpCnt)),col = "orange",lwd=2.5)

# mSmpT = table(cut(data$gintSmpCnt,breaks=c(c(0,5),seq(10,max(data$mutReg)+1, by = 10) )))
# mSmpT = rbind(mSmpT, round(mSmpT/sum(mSmpT) * 100,2))
# mSmpT = rbind(mSmpT, cumsum(mSmpT[2,]))
# rownames(mSmpT) = c("Freq","Percent","CumsumPercent")
# grid.newpage()
# grid.table(mSmpT, show.rownames = TRUE)


# # hist(data$mutReg, freq=F,breaks=20,font = 2,
#      xlab = "Number of driver regulators", ylab = "Frequency", main = "Distribution of Mutated Regulators",
#      col = "lightblue", border="gray")
# lines(density(as.numeric(data$gintSmpCnt)),col = "orange",lwd=2.5)

pdf(paste(figD,"barplot_mutatedSigRegulator_Apri21-2014.pdf",sep=""))

mRegT = table(cut(data$mutReg,breaks=c(c(0,5),seq(10,80,by=10),c(90,max(data$mutReg)+1) )))
mRegT = rbind(mRegT, round(mRegT/sum(mRegT) * 100,2))
mRegT = rbind(mRegT,cumsum(mRegT[2,]))
rownames(mRegT) = c("Freq","Percent", "CumSumPercent")
grid.newpage()
grid.table(mRegT, show.rownames = TRUE)
par(mar=c(6,6,2,1),mgp=c(4.5,0.6,0.5))
barplot(mRegT[2,], las=2,
        xlab = "Regulator Number", ylab = "%", main = "Regulators Number Distribution", font = 2,
        col = "lightblue", border="gray")
text(seq_len(length(mRegT[2,])) * 1.2 -  0.5 ,vapply(mRegT[2,], FUN=function(x){max(5,x*0.5)}, 1),
     col="black",font=2,cex=1.5,
     labels = mRegT[1,])


dev.off()

# 
# uniNumReg = unique(data$mutReg)
# uniNumReg = uniNumReg[order(uniNumReg)]
# alpha = 0.8
# plot(uniNumReg, vapply(uniNumReg,FUN=function(x){log(choose(x,ceiling(alpha * x)))}, 1),type = "l")
# plot(uniNumReg, vapply(uniNumReg,FUN=function(x){(choose(x,ceiling(alpha * x)))}, 1))
# 
# 
# uniNumSmp = data$gintSmpCnt
# cumsum(table(uniNumSmp))/sum(table(uniNumSmp))
# uniNumSmp = uniNumSmp[order(data$gintSmpCnt)]
# alpha = 0.8
# plot(uniNumSmp, vapply(uniNumSmp,FUN=function(x){log(choose(x,ceiling(alpha * x)))}, 1),type = "l")
# plot(uniNumSmp, vapply(uniNumSmp,FUN=function(x){min(500000,(choose(x,ceiling(alpha * x))))
#                                                  }, 1), pch = 19, col ="blue")
# abline(v=30,col="red")
# 
# plotCombPropSmp = function()
#   ## cover 80% of the 
#   temp = cumsum(table(uniNumSmp))/sum(table(uniNumSmp))
#   which.min(temp[which(temp > 0.8)])
