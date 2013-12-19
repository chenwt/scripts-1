
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


plotBigTable = function(gtInfo,rawcount,annoInfo,calledInfo,scaleCell,plotStep,title){
  ##Desp: plot mutations data 
  #Display: final mutation: 0/1 binary in calledInfo dataset
  #       raw maf infor: alt:ref in rawcount dataset(only output those in mutation cells)
  #       novel mutations in relapse?       
  plotStep = ifelse(plotStep =="", 0, as.numeric(plotStep))

  # table size
  numCol = ncol(rawcount)
  numRow = nrow(rawcount)
  
  # background canvas
  blankPlot = function(cellScale){  
    plot(c(-3*cellScale,-3*cellScale),col="white", 
         ylim=c(-3*cellScale,cellScale*(numRow+10)),
         xlim=c(-3*cellScale,cellScale*(numCol+10)),
         xaxt="n",yaxt="n",type="n",frame.plot=FALSE,
         xlab="",ylab="", xaxs="i",yaxs="i",mar=c(0, 0, 0, 0))
  }

  #data
  cellBgColors = c("white","orange")

  ##label and header
  patientType = sapply(colnames(rawcount),function(x){unlist(strsplit(x,"_"))[2] })
  patientID   = sapply(colnames(rawcount),function(x){unlist(strsplit(x,"_"))[1] })
  for (i in seq(1,length(patientID),by=3)) patientID[c(i+1,i+2)] = ""
  chrPos   = rownames(gtInfo)
  annoFuncColors = c("mediumpurple3","maroon1")
  annoFuncPlot   = ifelse(as.vector(sapply(annoInfo$FunctionClass,as.character)) == "nonsynonymousSNV",annoFuncColors[1],annoFuncColors[2])

  #text maf infor
  if(plotStep > 1) {
    #maf
    getAltAllele = function(x){as.numeric(unlist(strsplit(x,"/"))[1])}
    cellTextColors = c("white","blue")

    getPlotMaf1 = function(rawcount,calledInfo){
      numCol=ncol(rawcount);numRow=nrow(rawcount)
      out = matrix("",ncol=numCol,nrow=numRow)
      for (i in seq(1,numRow)){
        for (j in seq(1,numCol,by=3)){
          #print(calledInfo[i,(j+2)])
          if(calledInfo[i,(j+1)] > 0 || calledInfo[i,(j+2)] > 0){
            out[i,j:(j+2)] = sapply(rawcount[i,j:(j+2)],function(x){
              gsub("/",":",as.character(x))})
          }    
        }
      }
      return(out)
    }    
    rawcountPlotMaf1 = getPlotMaf1(rawcount,calledInfo)
    
    getPlotMaf2 = function(rawcount,calledInfo){
      totalAlt = apply(rawcount,1,function(x){as.numeric(unlist(strsplit(x,"/"))[1])})
      out = ifelse(totalAlt < 1,"",x)
      return(out)
    }
  }
    #rawcountPlotMaf2 =getPlotMaf2(rawcount,calledInfo)
    
    # rawcountPlotPene = as.data.frame(matrix(NA,nrow=numRow, ncol=numCol))
    # mafBorderColors = c('gray','red','green')
    # rawcountPlotPene = t(apply(rawcount,1,function(x){
    #   out = x
    #   for (i in seq(1,numCol,by=3)){
    #     n = getAltAllele(x[i+2]); t= getAltAllele(x[i+1]); r = getAltAllele(x[i+2])
    #     if ( t > 2 && r > 2){
    #       out[i] = out[i + 1] = out[i+2]=mafBorderColors[2]
          
    #     }else if (t <= 2 && t >2){
    #       out[i] = out[i+1] = out[i + 2] = mafBorderColors[3]
    #     }else{
    #       out[i] = out[i+1] = out[i+2] =  mafBorderColors[1]   
    #     }
    #   }
    #   return(out)}))


   #genotype
  if(plotStep > 2){
      gtValType  = as.character(unlist(unique(unlist(gtInfo))))
      gtColors   = c("white","yellow","orange")
      gtInfoColor = apply(gtInfo,c(1,2),as.character)
      for (i in 1:length(gtValType)){ gtInfoColor = gsub(gtValType[i], gtColors[i] ,gtInfoColor)} 
    }


  cellTextColors=c("white","blue")
  calledInfoTxt =ifelse(calledInfo==cellBgColors[1],cellTextColors[1],cellTextColors[2])
  
  pdf(paste(title, plotStep,".pdf",sep=""),
      height=numRow/numCol*40,width=40, onefile=TRUE, 
      family='Helvetica', pointsize=12)
  scaleCell = 10
  blankPlot(scaleCell)
  #plotting 
  for (ix in  scaleCell *(1:numCol)){
    for (iy in scaleCell *(1:numRow)){
      #mutation or not: fill color
      calledInfoTemp = calledInfo[(numRow - iy/scaleCell + 1),ix/scaleCell]
      rect(ix, iy, scaleCell + ix,scaleCell +iy,
           border="gray",lwd=0.5,
           col=ifelse(calledInfoTemp == 0,cellBgColors[1],cellBgColors[2]))          
      if(plotStep > 1){
        #maf info: text and color
        labelsTemp = rawcountPlotMaf1[(numRow-iy/scaleCell +1),ix/scaleCell]
        text(y=mean(iy,scaleCell + iy)+ 1/2 * scaleCell,x=mean(ix,scaleCell+ix)+ 1/2 * scaleCell,
             #labels = rawcountPlotMaf1[iy/scaleCell,ix/scaleCell],col=mafFillColors[1],
             labels = labelsTemp,
             #col   = ifelse(labelsTemp == "", cellTextColors[1],cellTextColors[2]),
             font=2) #,cex = 1.1
      }
      #if(plotStep > 2){
        
        # text(y=mean(iy,scaleCell + iy)+ 1/2 * scaleCell,x=mean(ix,scaleCell+ix)+1/2 * scaleCell,
          # labels = rawcountPlotMaf2[iy/scaleCell,ix/scaleCell], #col=mafFillColors[2], 
           # col=calledInfo[iy/scaleCell,ix/scaleCell],font=2) #,cex = 1.1
        #mutation penetration: border color
      #  rect(ix,iy,scaleCell + ix, scaleCell + iy,
       #     border = rawcountPlotPene[iy/scaleCell,ix/scaleCell], lwd = 5)
      #}
    }
  }
  
  #adding axies information
  #column header
    text(y=scaleCell* (numRow + 1.5) ,x=scaleCell*(1:numCol)+0.5 *scaleCell,
         labels=patientType,pos=3, font=2,,cex=1.5) #
         text(y=scaleCell* (numRow + 3 ) , x=scaleCell*(1:numCol) + 1.5*scaleCell,
         labels=patientID,pos=3,font=2,cex=1.5) #
    #right row label
    text(x=scaleCell * (numCol + 2) , y = scaleCell*((numRow+0.5):(1+0.5)) , 
        labels=as.vector(sapply(annoInfo$Gene,as.character)), 
        col   = annoFuncPlot,
        pos=4, cex=1.5, font=2) #
    text(x=scaleCell *(numCol + 5), y = scaleCell*((numRow+0.5):(1+0.5)) ,
        labels=chrPos,font=2, cex=1.5,
        col   = annoFuncPlot,
        ) #,
       
  #seperate patient by vertical line
    abline(v=scaleCell*seq(1,numCol+1,by=3),lwd=2)
  
  # #add legend
  #   legTxt = c(c("Not called","Called mutation"),"",
  #            paste("raw MAF:",c("no","called")), "",
  #            paste("mutType: ",c("No","Both Tumor&Relapse","Novel in Relapse")),"",
  #            paste("funcClass:",c("nonsynonymousSNV","stopgainSNV")))
  #   legColors = c(cellBgColors, "white", 
  #               mafFillColors, "white", 
  #               mafBorderColors, "white",
  #               annoFuncColors)
  #   ltyType = c(rep(NA,(length(gtColors) +1+ length(mafFillColors))),NA,
  #             rep(1,(length(mafBorderColors) +1+ length(annoFuncColors))))
  #   lwd = rep(2.5,length(ltyType))
  #   pchType = c(rep(15,(length(gtColors) +1+ length(mafFillColors))),NA,
  #             rep(NA,(length(mafBorderColors) +1+ length(annoFuncColors))))
  #   ptBgColors = legColors   
  # #legend(y = scaleCell * -1, x = scaleCell *  numCol / 2,
  # legend( y=scaleCell* (numRow + 8) , 
  #         x=scaleCell * (numCol + 2),
  #         #x = - 2 *scaleCell,
  #         text.font =2, box.col = "white", cex=1.5, #horiz = TRUE,
  #         fill = "white",border = "white",
  #         legend = legTxt,
  #         lty    = ltyType,       
  #         lwd    = lwd,
  #         pch = pchType,
  #         #pt.lwd = 5,
  #         pt.cex = 2.5,
  #         pt.bg = ptBgColors,
  #         col=legColors,
  # ) 
  box("figure",col="gray")
  dev.off()
} 

plotBigTable(gtInfo,rawcount,annoInfo,calledInfo,10,1, "imagePlotAllMutations_step")
plotBigTable(gtInfo,rawcount,annoInfo,calledInfo,10,2, "imagePlotAllMutations_step")

#plot scatter plot for one patients
rawcountPlot = apply(rawcount,c(1,2),function(x){
  temp = as.character(x)
  alt = as.numeric(unlist(strsplit(temp,"/"))[1])
  ref = as.numeric(unlist(strsplit(temp,"/"))[2])
  return(alt/ref)
  })

pdf("scatterPlotFinalMAF.pdf")
layout(matrix(1:16,byrow=T,nrow=4))
for (i in seq(1,48,by=3)){
  print(paste(i+1,i+2))
  plot(rawcountPlot[,c(i+1,i+2)] * 100, 
       pch="*",col="blue",cex=1.5,
       xlim=c(0,100),ylim=c(0,100),
       mar=c(0.1, 0.1, 0.1, 0.1))
}
dev.off()

##pos control
gene = read.table("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/reports/result_labmetgNov/fig/positiveControlGene.txt",sep="\t")
posGene = gene$V1
posGeneIndex = annoInfo$Gene %in% posGene == TRUE
rawcount = allData$rawcount[posGeneIndex,]
calledInfo = allData$calledInfo[posGeneIndex,]
gtInfo = allData$gtInfo[posGeneIndex,]
annoInfo = allData$annoInfo[posGeneIndex,]

plotBigTable(gtInfo,rawcount,annoInfo,calledInfo,10,2,"imagePlotPositiveControlGene_step")


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
