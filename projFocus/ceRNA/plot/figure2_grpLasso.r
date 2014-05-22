##after running regression to indentify candidate driver ceRNA, plot heatmap of target-regulator
## to show the predictibility of ceRNA regulators
###----under development----
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
figd = paste(rootd,"/DATA/projFocus/report/May2014/fig/",sep="")
cdt = paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="-")


##---load data
drRegfile = paste(rootd, "/DATA/projFocus/result/05012014/candiReg/runApr30/plotData/RYK_candidateRegs.txt",sep="")
expfile =  paste(rootd, "/DATA/projFocus/result/05012014/candiReg/runApr30/plotData/RYK_candidateRegs.txt_reg_exp_tumor.temp",sep="")

dataExp = read.table(expfile,sep="\t",header=T, stringsAsFactors=F)
allGenes = dataExp[,1]
dataExp = dataExp[,-1]; rownames(dataExp) = allGenes
tgene  = as.character(allGenes[1]); allRegs = allGenes[-1]

con  <- file(drRegfile, open = "r"); drRegsL <- NA 
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  templine = unlist(strsplit(oneLine,"\t"))
  if ( length(grep('^[[:digit:].e-]+$', templine[3])) > 0  & 
         as.numeric(as.character(templine[3]) ) <= 0.01) { 
    drRegsL = c(drRegsL,templine[1])
  } 
}
close(con); drRegsL = drRegsL[-1]

###--prepare data----
expRegD = apply(dataExp[c(drRegsL,tgene),],c(1,2),as.numeric)
expRegD = expRegD[,order(expRegD[tgene,])]
expRegD =t(apply(t(expRegD), 2, scale))
colnames(expRegD) = colnames(dataExp)
###-------------plot by heatmap.2----
require(gplots)
range(expRegD)

myCol = colorRampPalette(c("blue","white","red"))(n=299)
myBreaks = c(seq(-10,-1,length=100),seq(-1,1,length=100),seq(1,5.5,length=100))
lmat = rbind(c(0,3),c(2,1),c(0,4)) ; lwid = c(0.5,4.2); lhei = c(1,4,0.5)
par(cex.main = 2.5, )
heatmap.2(expRegD, Rowv=F,Colv=F,key=F,trace = "none",dendrogram="none",labCol="", font = 2,cex.main = 2.5,
          xlab= "Tumor samples",
          main = paste(tgene,"driver CeRNA dirvers from group lasso"),
          symm=F,symkey=F,symbreaks=T,breaks = myBreaks,
          lmat = lmat, lwid = lwid, lhei = lhei,
          scale='none',col = myCol)



###------plot by rect -----

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

###-----plot by ggplot2----
require(ggplot2)
require(reshape2)
require(grid)

datfx <- data.frame(indv=factor(paste("ID", 1:20, sep = ""),
                                levels =rev(paste("ID", 1:20, sep = ""))), matrix(sample(LETTERS[1:7],80, T), ncol = 4))
# converting data to long form for ggplot2 use
datf1x <- melt(datfx, id.var = 'indv')
plotx <-  ggplot(datf1x, aes(indv, variable)) + 
          geom_tile(aes(fill = value),colour = "white")  +   
          scale_fill_manual(values= terrain.colors(7))+ scale_x_discrete(expand=c(0,0))
px <- plotx

#Y axis quantitaive ggplot data
datfy <- data.frame(indv=factor(paste("ID", 21:40, sep = ""),
                                levels =rev(paste("ID",21:40, sep = ""))), matrix(sample(LETTERS[7:10],100, T), ncol = 5))
# converting data to long form for ggplot2 use
datf1y <- melt(datfy, id.var = 'indv')
ploty <-  ggplot(datf1y, aes( variable, indv)) + geom_tile(aes(fill = value),
                                                           colour = "white")  +   scale_fill_manual(values= c("cyan4", "midnightblue", "green2", "lightgreen")) + scale_x_discrete(expand=c(0,0))
py <- ploty  +  theme(legend.position="left",  axis.title=element_blank())

# plot XY quantative fill
allRegs = row.names(expRegD)
allSmps = colnames(expRegD)
datfxy <- data.frame(indv=factor(paste("ID", 1:20, sep = ""),
                                levels =rev(paste("ID", 1:20, sep = ""))), 
                     matrix(rnorm (400, 50, 10), ncol = 20))
datfxy = data.frame(gene=factor(allRegs,levels=allRegs),expRegD)
names (datfxy) <- c("gene",allSmps)
datfxy <- melt(datfxy, id.var = 'gene')
levels (datfxy$ variable) <- rev(allSmps)
pxy <- plotxy <-  ggplot(datfxy, aes(gene, variable)) + 
  geom_tile(aes(fill = value),colour = "white")  + 
  scale_fill_gradient2(low="blue",mid= "white",high="red",
                       breaks=c(seq(-10,-1,length=100),seq(-1,1,length=100),seq(1,5.5,length=100)),
                       guide = "colorbar") + 
  theme(axis.title=element_blank())

#Define layout for the plots (2 rows, 2 columns)
layt<-grid.layout(nrow=2,ncol=2,heights=c(7/8,1/8),widths=c(2/8,6/8),default.units=c('null','null'))
#View the layout of plots
grid.show.layout(layt)

#Draw plots one by one in their positions
grid.newpage()
pushViewport(viewport(layout=layt))
print(py,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(pxy,vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(px,vp=viewport(layout.pos.row=2,layout.pos.col=2))

