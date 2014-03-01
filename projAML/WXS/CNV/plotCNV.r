##~/tools/R/R_current/bin/Rscript
#-----function
chr2Num = function(chr) {
  if (chr == 'X' || chr =='x') 
  {chr = 23}
  else if (chr == 'Y' || chr == 'y')
  { chr = 24 }
  else if (length(grep(chr,vapply(1:22,as.character,'1'))) > 0 ) 
  { chr = as.integer(chr)}
  else 
  {print("error ")}
  return(chr)
}
blankPlot = function( maxX, maxY){  
  plot(c(-10,0),col="white", 
       xlim = c(0, maxX*1.05), ylim = c(0, maxY+5),
       xaxt="n",yaxt="n",type="n",frame.plot=FALSE,
       xlab="",ylab="", xaxs="i",yaxs="i",mar=c(0, 0, 0, 0))
}
cnv2color = function(str){
  if (str == "DEL") { str = "blue"}
  else if (str == "DUP") {str = "red"}
  else {print("error")}
  return(str)
}
##----endFunc

cwd = "/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/CNV/report/"
setwd(cwd)
figDir="/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/reports/result_Feb2014/fig"
cnv = paste(cwd,"/combined.xcnv_sorted",sep="")

chrInoffile = "/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/database/humanChromInfo.tsv"
##chrom information
chromLenInfo = read.table(chrInoffile,sep="\t",header=T,stringsAsFactors=F)
chromLenInfo = chromLenInfo[-nrow(chromLenInfo),]
for ( i in 1:nrow(chromLenInfo)){
  chromLenInfo$chr[i] = chr2Num(chromLenInfo$Chromosome[i])
  chromLenInfo$len[i]  = as.numeric(gsub(",","",chromLenInfo$Base.pairs[i]))  
}
coory = cumsum(chromLenInfo$len) 
names(coory) = chromLenInfo$Chromosome

##--
data = read.table(cnv,header=T,stringsAsFactors=F)
pts = unique(data$SAMPLE)
pts = pts[-c(grep("14",pts),grep("10",pts))]
numPts = length(pts)
numChr = 24

pdf(paste(figDir,"/plotCNV_02212014.pdf",sep=""),width=12)
par(mar=c(1,8,0.5,0))
blankPlot(max(coory),numPts)
##get unique pateints 
ptsUni = unique(vapply(pts,FUN=function(x){substr(x,11,16)},'a'))
##--highlight each pt bg
for (i in 1:length(pts)){
  temp_pid = substr(pts[i],11,16)
  colBg = ifelse(grep(temp_pid,ptsUni) %% 2 , "lightblue","white")
  rect(0, i-1, max(coory)+10000, i ,border="white",col=colBg)    
}
##--sep each chrom
for (i in 1:(length(coory))){ 
  if (i == 1) {
    lines(x=c(coory[i],coory[i]), y = c(0,numPts))
    text(coory[1]/2, numPts + 0.5, labels=names(coory)[1],cex=0.8)
  } else {
    lines(x=c(coory[i],coory[i]), y = c(0,numPts))
    text((coory[i-1] + coory[(i)] )/2, numPts + 0.5, labels=names(coory)[i],cex=0.8)
  }
}
##---axis label
axis(2,at=0.5:(length(pts)-0.5),cex.axis=0.7,labels=vapply(pts,FUN=function(x){
    substr(x,11,19)},'a'),las=2,lwd=0.3, tck = -0.01,  col.ticks="gray")
axis(1, at= max(coory)/2,labels="Replase: 04 Diagnosis: 03")

##--plot CNV data
for ( j in 1:nrow(data)){      
      temp_y = grep(data$SAMPLE[j],pts)
      if (length(temp_y) > 0 ){
            temp = unlist(strsplit(data$INTERVAL[j],"[:|-]"))
            temp_chr = chr2Num(temp[1])           
            if (temp_chr == '1') { 
              temp_xs = as.integer(temp[2] )
              temp_xe = as.integer(temp[3] )
            }else{
              temp_xs = as.integer(temp[2] ) + coory[temp_chr-1]
              temp_xe = as.integer(temp[3] ) + coory[temp_chr-1]
            }
            temp_col = cnv2color(data$CNV[j])
            rect(temp_xs+10000, temp_y-1, temp_xe+10000, temp_y,border=temp_col,col=temp_col)    
      }else{
            next
      }
}
dev.off()

getwd()
