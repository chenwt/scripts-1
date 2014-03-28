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
cdt = paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="-")

tf = paste(rootd, "DATA/projFocus/result/03102014/som/brca_somTumorWU.mat_gint_2k_dist0_Mar-11-2014.matrix",sep="")
tm = read.table(tf,header=T)
moccur = apply(tm[,-c(1:4)],1,function(x){length(x[which(x!=0)])})
data = table(cut(moccur,breaks=c(0,1,2,3,5,6,10,15,20,30,73)))/nrow(tm)*100

pdf(paste(figd,"/barplot_PromoterRegion_freqs_",cdt,".pdf",sep=""))
par(mfrow=c(1,1),mar=c(5,5,4,2),mgp=c(3,0.5,0))
barplot(data, las = 2, width=1,space=0.2,
        col = "lightblue",border="gray",
        ylab = "percentage (total#33271)",
        xlab="mutation recurrence",
        main = "Promotor region mutation's freq \n(Gint target genes, 73 samples)")
text(seq(from=(0.2+0.5),length.out=10,by=1.2), rescale(data,to=c(10,60)),
     labels=paste(round(data,2),"%",sep=""),cex=0.75,font=2)

dataR = sort(table(tm$geneName),decreasing=T)
data = table(cut(dataR,breaks=c(seq(0,500,by=50))))/nrow(dataR) * 100
barplot(data,
# barplot(data, 
     col="lightblue", las=2,border="gray",
      ylab= "percent of gene", width=1,space=0.2,
     main = " promotor region (+-2kb) mutations per gene")
text(seq(from=(0.2+0.5),length.out=10,by=1.2), rescale(data,to=c(5,25)),
     labels=paste(round(data,1),"%",sep=""),cex=0.75,font=2)
dev.off()



####-----freq plot for KIF23

ff = paste(rootd,"/DATA/projFocus/result/02022014/som/test/KIF23_mut_dist0b.matrix",sep="")
ddf = read.table(ff)
tempidx = vapply(smps,FUN= function(x){grep(x,colnames(rawSom))},1)
dataf = apply(ddf[,tempidx],c(1,2),as.numeric)
dataf = dataf[apply(dataf != 0,1,any), , drop=FALSE]
freq1b = apply(dataf,1,function(x){length(x[x!=0])})

pdf(paste(figd,"barplot_KIF23_mutRecurence",cdt,".pdf",sep=""))
par(mgp=c(2,1,0),mar=c(3.5,3.5,4,0))
dp = table(freq1b)/length(freq1b) * 100
barplot(dp,
        ylim = c(0,100),
        col="lightblue",border="gray",
        width=1,space=0.2,font=2,
        xlab = "point mutation recurrence", 
        ylab = paste("percentage < total#",length(freq1b),">",sep=""),
        main = "KIF23 point mutation recurence \n <42 samples> ")
# text(seq(from=(0.2+0.5),length.out=length(dp),by=1.2), rescale(dp,to=c(20,60)),
#      labels=round(dp,1),cex=0.75,font=2)      

###dignostic plot
#freqs 1kb
##dataSom from reportMar13_Plot_2.R
freq1k = apply(dataSom,2,function(x){length(x[x!=0])})
par(mgp=c(2,1,0),mar=c(3.5,3.5,4,0))
dp = table(freq1k)/length(freq1k) * 100
barplot(dp, ylim=c(0,40),
        col="lightblue",border="gray",
        width=1,space=0.2,
        xlab = "mutation hotspot recurrence", 
        ylab =  paste("percentage < total#",length(freq1k),">",sep=""),
        main = "KIF23 mutation hotspot recurence \n <42 samples> ")
# text(seq(from=(0.2+0.5),length.out=length(dp),by=1.2), rescale(dp,to=c(20,60)),
#      labels=round(dp,1),cex=0.75,font=2)      

dev.off()