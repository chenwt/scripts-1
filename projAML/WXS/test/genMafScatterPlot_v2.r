
sysInfo = Sys.info()
ifelse(sysInfo['sysname'] == "Darwin",
  setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/result_anno_filter_final/somFinal/"),
  setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/result_anno_filter_final/somFinal/")
)

#TODO:
# make three pdfs, show, mutations; maf of pateints with mutation; example of positive control genes
#raw count data
getData = function(){
  rawcount = read.table("All_patients_AlleleFreq_rawcount_sorted2.txt",sep="\t")
  colnames(rawcount) = sapply(rawcount[1,],as.character)
  rownames(rawcount) = rawcount[,1]
  rawcount = rawcount[,-1]
  rawcount = rawcount[-1,]  

  #genotype data
  gtInfo = read.table("All_patients_AlleleFreq_gt_sort2.txt",sep="\t")
  colnames(gtInfo) = sapply(gtInfo[1,],as.character)
  rownames(gtInfo) = gtInfo[,1]
  gtInfo = gtInfo[,-1]
  gtInfo = gtInfo[-1,]

  #annotation data
  annoInfo = read.table("allPatientsMutationAnnot.txt",sep="\t",header=T)

  #calledMut info
  calledInfo = read.table("allPatientsMutCalledSorted2.txt",sep="\t",header=T)
  rownames(calledInfo) = calledInfo[,1]
  calledInfo = calledInfo[,-1]
  return(list(rawcount=rawcount,
              gtInfo=gtInfo,
              annoInfo=annoInfo,
              calledInfo=calledInfo))
}
##prepare plot dataset
allData = getData()
rawcount = allData$rawcount#[1:20,]
calledInfo = allData$calledInfo#[1:20,]
gtInfo = allData$gtInfo#[1:20,]
annoInfo = allData$annoInfo#[1:20,]


#plot scatter plot for one patients
rawcountPlot = apply(rawcount,c(1,2),function(x){
  temp = as.character(x)
  alt = as.numeric(unlist(strsplit(temp,"/"))[1])
  ref = as.numeric(unlist(strsplit(temp,"/"))[2])
  return(alt/(alt+ref))
  })

##---debuging PARVUA tumor relapse order
targetIndex = grep("PARVUA", colnames(rawcountPlot))
rawcountPlot = rawcountPlot[,c(1:37,39,38,40:48)]

##------plotting
pdf("scatterPlotFinalMAF.pdf")
layout(matrix(1:16,byrow=T,nrow=4))
par(mar=c(2.8,2.8,0.5,0.5),mgp=c(1.5,0.5,0))
for (i in seq(1,48,by=3)){
  print(paste(i+2,i+1))
  plot(rawcountPlot[,c(i+2,i+1)] * 100, 
       pch="*",col="blue",cex=1.5,
       xlim=c(0,100),ylim=c(0,100)
       )
}
dev.off()

#----plotting with color---

rawcountPlotTumor = apply(rawcountPlot,c(1,2),function(x){
    temp = as.character(x)
    alt = as.numeric(unlist(strsplit(temp,"/"))[1])
    ref = as.numeric(unlist(strsplit(temp,"/"))[2])
    return(alt/ref)
  })

rawcountPlotNormal = rawcountPlot[,seq(1,48,by=3)] * 100
rawcountPlotTumor  = rawcountPlot[,seq(2,49,by=3)] * 100
rawcountPlotRelapse  = rawcountPlot[,seq(3,50,by=3)] * 100
patientID = sapply(colnames(rawcount)[seq(1,48,by=3)],function(x){substr(x,0,6)})

pdf("scatterPlotFinalMAFColored.pdf")
layout(matrix(1:16,byrow=T,nrow=4))
par(mar=c(2.8,2.8,0.5,0.5),mgp=c(1.5,0.5,0))
for (i in seq(1,16,by=1)){
  print( patientID[i])
  colors = rep("blue",length(rawcountPlotRelapse[,i]))
  colors[intersect(which(rawcountPlotRelapse[,i] > 10),which(rawcountPlotTumor[,i] < 10))] = "red"
  colors[intersect(which(rawcountPlotRelapse[,i] < 10),which(rawcountPlotTumor[,i] > 10))] = "green"
  plot(x = rawcountPlotRelapse[,i],rawcountPlotTumor[,i], 
       pch="*",col=colors,cex=2.0,
       xlab = paste(patientID[i],"Relapse"), ylab=paste(patientID[i],"Tumor"),
       xlim=c(0,100),ylim=c(0,100)
  )
}
dev.off()

##---output relapse specific mutations
out = rep("",6)
chrPos = rownames(rawcountPlot)
for (i in seq(1,16,by=1)){
   pid = patientID[i]
   comment = "relapse specific"
   info = annoInfo[intersect(which(rawcountPlotRelapse[,i] > 10),which(rawcountPlotTumor[,i] < 10)),]
   if (dim(info)[1] > 0) {
        for (j in 1:nrow(info)) {out = rbind(out, c(pid,comment,sapply(info[j,],as.character)))}
   }
  #annoInfo[intersect(which(rawcountPlotRelapse[,i] < 10),which(rawcountPlotTumor[,i] > 10))]
}
colnames(out) = c("patientID","comments", colnames(annoInfo))
write.table(out,"tableSomaticMutationRelapseSpecific.txt",col.names=F,sep="\t",quote=F)


##pos control
gene = read.table("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/reports/result_labmetgNov/fig/positiveControlGene.txt",sep="\t")
posGene = gene$V1
posGeneIndex = annoInfo$Gene %in% posGene == TRUE
rawcount = allData$rawcount[posGeneIndex,]
calledInfo = allData$calledInfo[posGeneIndex,]
gtInfo = allData$gtInfo[posGeneIndex,]
annoInfo = allData$annoInfo[posGeneIndex,]

plotBigTable(gtInfo,rawcount,annoInfo,calledInfo,10,2,"imagePlotPositiveControlGene_step")

###genes with high penetrality
data = calledInfo[which(rowSums(calledInfo) >1),]
apply(data,1,function(x){names(x>0)})
for ( i in 1:nrow(data))
print (paste(i, unique(sapply(colnames(data)[which(data[i,] !=0)],function(x){substr(x,1,6)}))) )
annoInfo[annoInfo[,1] %in% rownames(data[c(4,19,24),]),]

#other function
#val2col<-function(z,col1="blue",col2="white",col3="red",nbreaks=30,type="discrete"){  
#   if(type == "continious"){
#     if(is.character(z)){
#       z<-as.numeric(as.factor(z))
#     }
#     extreme=round(max(abs(z)))
#     breaks <- seq(-extreme, extreme, length = nbreaks)
#     ncol <- length(breaks) - 1
#     col <- colorpanel(ncol,col1,col2,col3)
#     CUT <- cut(z, breaks=breaks)
#     colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
#     names(colorlevels)<-names(z)}
#   else if(type == "discrete"){
#     if(is.factor(unique(unlist(z)))){
#       
#     } 
#     
#   }
  
#  return(colorlevels)
#}
