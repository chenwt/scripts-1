##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep.11,2013 
#TODO: Plot the heatmap; do this to muTect Data
#Description:integrating analysis for 16 calling from samtools 

#### run the below twice: one for samtools and one for gatk
rootd <- "/ifs/scratch/c2b2/ac_lab/jh3283/"
# rootd <- "/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/"

##----------------------------
# setwd(paste(rootd,"projAML/WXS/callVars/report/",sep=""))
setwd(paste(rootd, "projAML/WXS/muTect/report/",sep=""))
fn <- list.files(pattern="RN-TN")
data <- lapply(fn, function(x) read.table(x,sep="\t",stringsAsFactors=F))
names(data) <- gsub("_RN-TN.txt","",fn)
# V2,V3, V22,v23,v24P

# samtools
chr_pos <- sapply(data,function(pt) paste(pt[,22],pt[,23],sep="_") )
names(chr_pos) <- names(data)
data_mat <- sapply(data,function(pt){
		  lapply(chr_pos,function(x){
		  tmp <-  unlist(strsplit(x,"_"))
		  pt[which(pt[,22]==tmp[1] & pt[,23]==tmp[2]),]}
      )})
		  
getFeatureList <- function(pt) {
  tmp <- paste(pt[,22],pt[,23],sep="_")
  names(tmp) <- names(data)
  return(tmp)	      }

list2mat <- function(inlist){ 
  nr <- unique(unlist(inlist))
  mat.temp <- lapply(inlist,function(x) nr %in% x * 1)
  mat <- data.frame(mat.temp)
  colnames(mat) <- names(inlist)
  rownames(mat) <- nr
  return(mat)
}

#get heatmap of mutations in different patients
cp_mat <- list2mat(chr_pos)
cp_mat$chr_pos <- rownames(cp_mat)

#get the pos_to_annotation matrix
for (i in 1:ncol(cp_mat)){
  x <- colnames(cp_mat)[i]
  tmp <- data[[x]]
  if (i==1) pos2gene <- data.frame(cbind(chr_pos=paste(tmp[,22],tmp[,23],sep="_"),gene=tmp[,2],func=tmp[,3]))
  else pos2gene <- rbind(pos2gene,data.frame(cbind(chr_pos=paste(tmp[,22],tmp[,23],sep="_"),gene=tmp[,2],func=tmp[,3])))
  }
pos2gene <- unique(pos2gene)


#exclude false positive genes
mat <- merge(pos2gene,cp_mat,by="chr_pos")

# FPgene <- as.character(unlist(read.table(paste(rootd,"ref/databaseSeq/seq_FalsePositiveGenes.txt",sep=""))))
# fpidx <- unlist(sapply(FPgene,function(x)grep(x,mat$gene)))

# ifelse(length(fpidx)>0, mat <- mat[-fpidx,],mat<-mat)

mutGene_twice <- unique(mat$gene[duplicated(mat$gene)])
plotMat <- mat[unlist(sapply(mutGene_twice,function(x)grep(x,mat$gene))),]
write.table(plotMat,"table_16patient_mutation.txt",row.names=F,sep="\t",quote=F)



# require(xlsx)
# # write.xlsx(Mat,file="table_16patient_mutation.xlsx")
# write.xlsx(plotMat,file="table_16patient_mutation.xlsx")

# source(paste(rootd,"scripts/myUseful/R/heatmap/heatmap3.R",sep=""))
# plotdata <- as.matrix(plotMat[,-c(1:3)])
# plotdata <- plotdata[order(rowSums(plotdata),decreasing=T),]
# labRow <- paste(plotMat$gene,plotMat$chr_pos)
# labCol <- colnames(plotMat[,-c(1:3)])

# heatmap(plotdata,col= cm.colors(256),labRow=labRow)

##----------------------------
#integrate samtools and gatk

mutg.temp <- read.table(paste(rootd,"projAML/WXS/muTect/report/table_16patient_mutation.txt",sep=""),
  header=T,sep="\t")
mutg.temp2 <- read.table(paste(rootd,"projAML/WXS/callVars/report/table_16patient_mutation.txt",sep=""),
  header=T,sep="\t")

mutg <- rbind(cbind(mutg.temp2,cmt=rep("gatk",nrow(mutg.temp2))),
              cbind(mutg.temp,cmt=rep("samtools",nrow(mutg.temp))))

##----------------------------
#Filtering with cosmic census


#filtering wiht TCGA NEJM, TCGA provisional

require(xlsx)
write.xlsx(mutg[order(mutg$chr_pos),],
  paste(rootd,"projAML/WXS/reports/table_16patient_mutation_twice_Gene.xlsx",sep=""))
write.table(mutg[order(mutg$chr_pos),],
  paste(rootd,"projAML/WXS/reports/table_16patient_mutation_twice_Gene.txt",sep=""),
  row.names=F,quote=F,sep="\t")

##----------------------------
#output chr pos pid to grap MAF information for Normal Tumor and Relapse
pid <- unlist(read.table(paste(rootd,"projAML/WXS/PID_16.txt",sep="")))
mutgChrPosPid <- cbind(mutg$chr_pos[which(mutg==1,arr.ind=T)[,1]],
     colnames(mutg)[which(mutg==1,arr.ind=T)[,2]])
write.table(mutgChrPosPid,paste(rootd,"projAML/WXS/reports/mutgChrPosPid.txt",sep=""),
    row.names=F,quote=F,col.names=F,sep="\t")


##----------------------------
#plot MAF changing for each gene

setwd("/Volumes/ac_lab/jh3283/SCRATCH/projAML/WXS/reports/result_labmetgSep24")
dataraw <- read.table("result_varFreq_v2_overlap_clean.txt",header=T)
# gene_name chr       pos    PID N_Total N_Alt T_Total T_Alt R_Total R_Alt
#  NRAS   1 115258744 PANTWV     396     0     226   104     295    83

g2plot <- unique(data$gene_name[duplicated(data$gene_name)])
plotd <- data[data$gene_name %in% g2plot,]
plotd <- data.frame(gene=plotd$gene_name,
      N=plotd$N_Alt/plotd$N_Total,T=plotd$T_Alt/plotd$T_Total,
      R=plotd$R_Alt/plotd$R_Total)
myplotG <- function(gene,d){
    # tempd <- subset(d,d$gene==gene), not working 
    tempd <- d[which(d$gene==gene),]
    # print(dim(tempd))
    pdf(paste("graph_varfreq_v2/",gene,".pdf",sep=""))
    plot(x=1:3,type="n", ylim=c(0,0.6),
      xlab="sample code", ylab="MAF",xaxt="n")
    axis(1, at= c(1,1.5,2,2.5,3),labels=c("Normal","","Tumor","","Relapse"))
    apply(tempd,1,function(x)lines(x=1:3,y=x[2:4],
        type="b",lwd=1.5,pch=1,col="blue"))
    title(paste(gene,"Minor Allele frequency"))
    dev.off()
}
myplotG(g2plot[1],plotd)
sapply(g2plot,function(x)myplotG(x,plotd))

##----------------------------
#plot all other genens in one plot

plotd <- data[!(data$gene_name %in% g2plot),]
plotd <- data.frame(gene=plotd$gene_name,
      N=plotd$N_Alt/plotd$N_Total,T=plotd$T_Alt/plotd$T_Total,
      R=plotd$R_Alt/plotd$R_Total)
ngenes <- nrow(plotd)


# Create Line Chart

# convert factor to numeric for convenience 
# plotd$gene <- as.numeric(plotd$gene) 

# get the range for the x and y axis 
myplotM <- function(x,plotd) {
  pdf(paste("graph_varfreq_v2/single_",x[1],".pdf",sep=""))
    plot(x=1:3,type="n", xlab="sample code", ylim=c(0,0.6),ylab="MAF" ) 
    lines(x=1:3,x[2:4],type="b",lwd=1.5,col="blue")
  dev.off()
} 

apply(plotd,1,function(x)myplotM(x,plotd))
# set up the plot 
plot(x=1:3,type="n", xlab="sample code", ylim=c(0,0.6),ylab="MAF" ) 
colors <- rainbow(ngenes) 
# plotchar <- seq(1,1+ngenes,1)

# add lines 
for (i in 1:ngenes) { 
  if(plotd[i,4] > plotd[i,3]) 
    lines(x=1:3,plotd[i,2:4],type="b",lwd=1.5,col=colors[i])
  else
    lines(x=1:3,plotd[i,2:4],type="b",lwd=1.5,col="gray")
    # i <- 1
    # apply(plotd,1,function(x){lines(x=1:3,y=x[2:4],type="b", lwd=1.5,
     # col=colors[i]); i = i+1 })   
} 
# add a title and subtitle 
title("MAF plot", "genes have 1 mytations")

# add a legend 
for (i in 1:ngenes) { 
  if(plotd[i,4] > plotd[i,3]) 
    legend(1:3, range(c(0,0,6)), 1:ngenes, cex=0.8, col=colors,
     pch=plotchar, lty=linetype, title="MAF plot")
  else
    legend(1:3, range(c(0,0,6)), 1:ngenes, cex=0.8, col=colors,
     pch=plotchar, lty=linetype, title="MAF plot")
}




##----------------------------
#integrate cnv and viper results
setwd("~/SCRATCH/projAML/WXS/reports/result_labmetgSep24/")
load("~ys2559/scripts/AML/viperPerSubclass.rData")
require(org.Hs.eg.db)
geneSym <- org.Hs.egSYMBOL
geneSym <- as.list(geneSym[mappedkeys(geneSym)])
rownames(amlVp) <- geneSym[rownames(amlVp)]
mrs <- unlist(read.table("temp_MRs4VP.txt",header=F))
pid <- tolower(unlist(read.table("../../PID_16.txt",header=F)))
idxr <- na.omit(unlist(sapply(mrs,function(x) match(x,rownames(amlVp)))))
idxc <- na.omit(unlist(sapply(pid,function(x) match(x,colnames(amlVp)))))
mrspidVP <- mrsVP[idxr,idxc]
# mrspidVP <- mrspidVP[,-ncol(mrspidVP)]


tcnv <- read.table("TN.xcnv.cnv.bed.geneAnnot.pidtgene",header=T,sep=" ")
colnames(tcnv) <- c("pid","type","gene")
require(reshape2)
tcnv.wd <- dcast(
            transform(tcnv,id=rep(1:length(unique(tcnv$pid)),times=table(tcnv$pid))),
            gene~id,value.var="type",fill=0)
colnames(tcnv.wd) <- tolower(c("gene",unique(tcnv$pid)))
rownames(tcnv.wd) <- tcnv.wd[,1]
tcnv.wd <- tcnv.wd[,-1]
idxc <-  na.omit(unlist(sapply(pid,function(y) match(y,colnames(tcnv.wd)))))
idxr <- na.omit(unlist(sapply(mrs,function(x) match(x,rownames(tcnv.wd)))))
mrstcnv.wd <- tcnv.wd[idxr,idxc]
mrstcnv.wd <- apply(mrstcnv.wd,c(1,2),function(x){ifelse(x >0, x<- x-2,x)} )
# mrstcnv.wd <- mrstcnv.wd[,-ncol(mrstcnv.wd)]


rcnv <- read.table("RN.xcnv.cnv.bed.geneAnnot.pidtgene",header=T,sep=" ")
colnames(rcnv) <- c("pid","type","gene")
require(reshape2)
rcnv.wd <- dcast(
            transform(rcnv,id=rep(1:length(unique(rcnv$pid)),times=table(rcnv$pid))),
            gene~id,value.var="type",fill=0)
colnames(rcnv.wd) <- tolower(c("gene",unique(rcnv$pid)))
rownames(rcnv.wd) <- rcnv.wd[,1]
rcnv.wd <- rcnv.wd[,-1]
idxr <- na.omit(unlist(sapply(mrs,function(x) match(x,rownames(rcnv.wd)))))
idxc <-  na.omit(unlist(sapply(pid,function(y) match(y,colnames(rcnv.wd)))))
mrsrcnv.wd <- rcnv.wd[idxr,idxc]

mrsrcnv.wd <- apply(mrsrcnv.wd,c(1,2),function(x){ifelse(x >0, x<- x-2,x)} )
# mrsrcnv.wd <- mrsrcnv.wd[,-ncol(mrsrcnv.wd)]

##----------------------------
#calculate patient level CNV_MRs relationship
idxc <- intersect(colnames(mrspidVP),intersect(colnames(mrsrcnv.wd),colnames(mrstcnv.wd)))
idxrt <- intersect(rownames(mrspidVP),rownames(mrstcnv.wd))
idxrr <- intersect(rownames(mrspidVP),rownames(mrsrcnv.wd))
mrsrcnv.wd[idxrr,idxc] * mrspidVP[idxrr,idxc]
mrstcnv.wd[idxrt,idxc] * mrspidVP[idxrt,idxc]








