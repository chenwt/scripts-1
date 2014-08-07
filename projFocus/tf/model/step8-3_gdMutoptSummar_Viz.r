
##initiate
rm(list=ls())
usage = "Usage: Rscript step3-4_greedyOptCorr.r  -exp <expression file from python> -mut <mutation file from python>"
ERR = "ERROR:"
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,4,7)],collapse="-")
timeStart = Sys.time()

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
figd = paste(rootd, "/DATA/projFocus/report/Jul2014/fig/",sep="")

source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

### test files
input_act = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/summary/all_targets_actSuumary"
input_perm = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/summary/all_targets_permSuumary"
output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/summary/all_targets_mutSampleVector"

dataAct = read.csv(input_act,sep="\t",stringsAsFactor = F)
rownames(dataAct) <- dataAct$X; dataAct=dataAct[,-1]; dataAct = data.frame(apply(dataAct,c(1,2),as.numeric))
dataPerm = read.csv(input_perm, sep = "\t", stringsAsFactor = F)
rownames(dataPerm) <-  dataPerm$X; dataPerm = dataPerm[,-1]; dataPerm = data.frame(apply(dataPerm, c(1,2), as.numeric))

##--func--
extrat_gPerm = function(g, permName){
  ###func to extrac indice of all permutations for one gene
  tempname = vapply(permName, FUN=function(x){unlist(strsplit(x,"_"))[1]},'a')
  return(which(tempname == g, arr.ind=T))
}

getPval_actVSperm = function(gAct, gPerm) {
  nperm = NROW(gPerm)
  
  before = max(length(which(gAct[,5] < gPerm[,5], arr.ind=T) )/nperm, 1/nperm)
  after1 = max(length(which(gAct[,7] < gPerm[,7], arr.ind=T) )/nperm, 1/nperm)
  after2 = max(length(which(gAct[,9] < gPerm[,9], arr.ind=T) )/nperm, 1/nperm)
  
  return(c(rownames(gAct), before, after1, after2))
}

###--func end--

gls = rownames(dataAct)

### pval_actVSperm_before

pval_actVSperm = data.frame(matrix(NA,nrow=NROW(dataAct),ncol=4)); colnames(pval_actVSperm) = c("gene", "before", "after1", "after2");
cnt = 0
con <- file(input_perm, open='r')
g_prev = '';gPerm = rep(NA, NCOL(dataAct)); htag = 0; stag = 0
while(length(line_crt <- readLines(con,n=1, warn=F))==1){
  if (htag == 0){ htag = 1 ;next  }
  if (stag == 0) {g_prev =unlist(strsplit(line_crt,"_"))[1]; stag = 1;next }
    
  g_crt = unlist(strsplit(line_crt,"_"))[1]
  line_crt = unlist(strsplit(line_crt, "\t")); name_crt = line_crt[1]; val_crt = line_crt[-1]
 
 if (g_crt == g_prev ) {
    gPerm = rbind(gPerm, val_crt )
    g_prev = g_crt
    
 }else{
   ## do function
   #print(c(g_prev,g_crt))
   cnt = cnt + 1 #; if(cnt > 20) break
   
   pval_actVSperm[cnt, ] = getPval_actVSperm(dataAct[g_crt,],gPerm)
   g_prev = g_crt
   gPerm = rep(NA, NCOL(dataAct))
   
 }
 
}
close(con)


  #     plot(x=3:(length(gPerm) + 2), y = -log10(gPerm), 
  #          main = g, xlab = "index", ylab = "-log10(pval_optVSact)",
  #          col="gray", pch=19)
  #     points(x = 2, y =-log10(gAct), col = "red", pch = 19)
  #     text(x=1, y=max(-log10(gPerm)), labels=paste("Pval(n=",nperm,"): ", pval_perm, sep=""),font=2,pos=4)

###output
write.table(pval_actVSperm[order(pval_actVSperm$before),], paste(output,".compareStats",sep=""), sep="\t",row.names = F, quote = F)



###---plot
pdata = data.frame(pval_actVSperm)
pdata2 = pdata[which(pdata$before<=pCut2 & pdata$after1<=pCut2),]  

pdf(paste(output,".pdf",sep=""))
  plot(x=-log10(dataAct$pval_act), y=-log10(dataAct$pval_actVopt))
  abline(h=2, col="red",lwd=2)
  
  plot(x=pval_actVSperm_before, y=pval_actVSperm_after,
       pch = 19, col ="darkgray",
#        xlim = c(1,max(-log10(pval_actVSperm_before)) + 2),
     xlab = "pval_before",ylab = "pval_after", main = "pval of act v.s perm ")
  points(pdata2,pch = 19, col="red")
  abline(v=0.01,col="red")
  abline(h=0.01,col="red")
dev.off()


##----analysis part
pCut1 = -log10(0.01)
pCut2 = 0.05
print(c("allgene","act_before"))
nrow(dataAct)
NROW(dataAct[which(-log10(dataAct$pval_act) >= pCut1),])
NROW(dataAct[which(-log10(dataAct$pval_opt) >= pCut1),])
NROW(dataAct[which(-log10(dataAct$pval_actVopt) >= pCut1),])


NROW(pval_actVSperm[which(pval_actVSperm$before<=pCut2),])
NROW(pval_actVSperm[which(pval_actVSperm$after1<=pCut2),])

NROW(pval_actVSperm[which(pval_actVSperm$before<=pCut2&pval_actVSperm$after1<=pCut2),])
