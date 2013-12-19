setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/reports/result_labmetgNov")
scatterPlot = function(X,Y,label,title){
  xyLimits = c(min(X,Y),max(X,Y))
  plot(X,Y,
       xlab = "Relapse patient", ylab ="Tumor patients", main = title,
       xlim = xyLimits, ylim = xyLimits,
       pch = 15, col = "blue"
  )
  abline(a=0,b=1,col="gray")
  text(x=X,y = Y *(1 - 0.02),labels=label,
       cex = 0.5)
}

data = read.table("fig/rawVariantCount.txt",sep="\t",header=T)
pdf("fig/scatterPlotRawMutation.pdf")
scatterPlot(Y=data$TuA,X=data$ReA,data$Patient,title="Raw Variannts Number")
dev.off()

data = read.table("fig/rawReadsNumbersCount.txt",sep="\t",header=T)
pdf("fig/scatterPlotrawReadsNumber.pdf")
scatterPlot(X=data$Relapse, Y=data$Tumor,label=data$PID,title="Raw Reads Number")
dev.off()

data =read.table("fig/finalVrariantCount.txt",sep="\t",header=T)
pdf("fig/scatterPlotFinalVrariantCount.pdf")
scatterPlot(X=data$Re,Y=data$Tu,label=paste(data$patient,data$Total,sep=":"))
debv.off()
barplot(as.matrix(data[,c(2,3)]),col=c("green","red"),
        #legend = rownames(data),
        #legend = data$patient,        
        beside=TRUE,
        border=NA,space=0.1, names.arg = data$patient)
