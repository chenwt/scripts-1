##after running regression to indentify candidate driver ceRNA, plot heatmap of target-regulator
## to show the predictibility of ceRNA regulators

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
figd = paste(rootd,"/DATA/projFocus/report/topDown_02042014/fig/",sep="")
cdt = paste(unlist(strsplit(date()," "))[c(2,4,6)],collapse="-")
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

###-------------plot
dataExpPlot = matrix(val2col(unlist(dataExp),col=myCol),nrow=nrow(dataExp), byrow=F)
rownames(dataExpPlot) = rownames(dataExp)
colnames(dataExpPlot) = colnames(dataExp)

drawVectorHeat = function(z, pos = 0, size = 0.5) {
  myCol = colorRampPalette(c("blue","white","red"))(256)
  z = val2col(z,col=myCol)
  for (i in 1:length(z)){
    rect(i-1,ybottom=pos ,xright=i,ytop=(pos+ size),col=z[i],border=z[i])  
  }
}

addPIDLabel = function(z){
  axis(1,at=seq(0.5,length(z),1),labels=z,cex.axis=0.8,font=2,
       las=2, lwd=0.3, tck = -0.01, col.ticks="gray",)
}

addGeneLable = function(label,pos){
  axis(2,at=pos,cex.axis=1,font = 2 ,
       labels = label,las=2,lwd=0.3, tck = -0.01,  col.ticks="gray")
}

##recorder....
dataExp = dataExp[order(dataExp[,1]),]
allgenes = colnames(dataExp)
# dataExp = dataExp[,c(allgenes[1],cands,allgenes[-1][is.na(match(allgenes[-1],cands))])]
dataExpSort = dataExp[,c(allgenes[1],cddts,allgenes[-1][is.na(match(allgenes[-1],cddts))]) ] ##group from grouping using kmeans

dataExpNorm = apply(dataExp,2,function(x){rescale(x,to=range(dataExp[,1]))})

pdf(paste(figd,"heatmap_CHEK1_candidateRegulator",cdt,".pdf",sep=""))
require("scales")
plotHeatmap <- function (data) {
  
  par(mar=c(6,6,1,0),mfrow=c(1,1))
  blankPlot(maxX=nrow(data),maxY=ncol(data))
  drawVectorHeat(dataExp[,1],pos=0,size=2)
  addGeneLable(pos=1,label=colnames(data)[1])
  addPIDLabel(rownames(data))
  for (i in 1:(ncol(data)-1)){
    temp = rescale(data[,(i+1)],to=range(data[,1]))
    templab = colnames(data)[i+1]
    size = 1
    drawVectorHeat(unlist(temp),pos=(i * size) + 2,size=size)
    addGeneLable(pos= i*size - size/2 + 3,label=templab)
  }
  title(main=paste("\n\n",colnames(data)[1]," Regulators"))
}
plotHeatmap(dataExpSort)

require(gplots)
mycol=colorRampPalette(c("blue","white","red"))(20)
heatmap.2(as.matrix(t(dataExpSort)),col=mycol,scale="none",trace="none",Rowv="none",Colv="none",dendrogram="none")

dev.off()
