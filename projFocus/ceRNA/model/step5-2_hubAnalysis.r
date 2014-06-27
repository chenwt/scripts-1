##given the number 

CWD ="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/"
figd = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jun2014/fig/"
cdt = paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="-")
setwd(CWD)

###---func---
require(igraph)
loadNetFromKeyregsumfile <- function (netfile) {
  ##load 2 column network file from keyRegsummary file
  dataDF = as.data.frame(matrix(NA, ncol=2))
  colnames(dataDF) = c("target", "regulator")
  con <- file(netfile, open='r')
  while (length(linecrt <- readLines(con, n =1, warn = FALSE)) > 0){
    vectcrt = unlist(strsplit(linecrt,"\t"))
    regCrt = unlist(strsplit(vectcrt[2],";"))
    dataDF <- rbind(dataDF, cbind(target = rep(vectcrt[1],length(regCrt)), regulator = regCrt))
  }
  close(con)
  return(unique(na.omit(dataDF)))
}
loadMutGene <- function (mutTFbsfile) {
  mutTF = vector(mode="character")
  con = file(mutTFbsfile, open="r")
  while (length(linecrt <- readLines(con, n = 1, warn = F)) > 0){
    mutTF <- c(mutTF, unlist(strsplit(linecrt, "\t"))[1])
  }
  close(con)
  mutTF <- mutTF[-1]
  return(as.data.frame(table(mutTF)))
}

##load network file
netfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_flex_max_1000.tsv.significant.summary.netfile"
netRaw = loadNetFromKeyregsumfile(netfile)
netRaw = netRaw[-1,]
allDrRegs = unique(netRaw[,2])


###
data2 = graph.data.frame(netRaw,directed=FALSE)
g = simplify(data2)
orderedVertex = get.data.frame(g,what = "vertices")
g = simplify(data2)

###---plot scatter plot
pdf(paste(figd,"/step5-2_hubAnalysis_allSigceRNADriver_", cdt,".pdf",sep=""))

allDrSort = sort(table(netRaw[,2][-1]),decreasing=T)
plot(0,xlim=c(0,200),ylim=c(0,100),type="n",
     xlab="Top ceRNA Driver \n(Rank by number of targets)",
     ylab="Percentage of total targets",
     font =2)

for (i in seq(10,200,by=10)){
  topDrRegs <- row.names(as.matrix(allDrSort[1:i]))  
  topDrNet = netRaw[netRaw[,2] %in% topDrRegs, ]
  
  #points(x=i, y=length(unique(topDrNet[,1]))/1086*100,col="red",pch=19,lwd=3,)
  print(paste(length(unique(topDrNet[,2])),length(unique(topDrNet[,1]))/1086*100))
}
dev.off()
write.table(file=paste(CWD,"/brca_greedy_sig_ceRNADriver_rankByNumOfTarget_",cdt,sep=""), as.matrix(allDrSort),sep="\t",quote=FALSE,col.names=F)


####----plot top 20 ceRNA network
i = 20
pdf(paste(figd,"/step5-2_hubAnalysis_hubceRNADriver_top",i,"net.pdf",sep=""))
topDrRegs <- row.names(as.matrix(allDrSort[1:i]))  
topDrNet = netRaw[netRaw[,2] %in% topDrRegs, ]

###--plot network of hub regulator
data2 <- graph.data.frame(as.data.frame(topDrNet),directed=FALSE)
g <- simplify(data2)
orderedVertex = get.data.frame(g,what = "vertices")
allVertexName <- orderedVertex$name
g <- simplify(data2)

# str(V(g))
# for (i in 1: length(allVertexName)){
#   pattern = paste("^", allVertexName[i], "$", sep="")
#   if (length(grep(pattern,mutTFgene)) > 0 & length(grep(pattern, mutMirgene)) > 0 ){
#     V(g)$color[i] <- "red"
#     V(g)$label[i] <- allVertexName[i]
#   }else if(length(grep(pattern,mutTFgene)) == 0 & length(grep(pattern, mutMirgene)) > 0){
#     
#     V(g)$color[i] <- "orange"
#     V(g)$label[i] <- allVertexName[i]
#     
#   }else if(length(grep(pattern,mutTFgene)) > 0 & length(grep(pattern,mutMirgene)) == 0) {
#     
#     V(g)$color[i] <- "steelblue"
#     V(g)$label[i] <- allVertexName[i]
#     
#   }else{
#     V(g)$color[i] = "gray"
#     V(g)$label[i] = ""
#     
#   }
# }

V(g)$color <- ifelse(vapply(orderedVertex$name,FUN=function(x){pattern = paste("^", x, "$", sep=""); length(grep(pattern,topDrRegs))},1) > 0,
                     "red", "lightblue")
V(g)$size <- ifelse(vapply(orderedVertex$name,FUN=function(x){pattern = paste("^", x, "$", sep="");length(grep(pattern,topDrRegs))},1) > 0,
                    2, 0.8)
V(g)$label.cex <- ifelse(vapply(orderedVertex$name,FUN=function(x){pattern = paste("^", x, "$", sep="");length(grep(pattern,topDrRegs))},1) > 0,
                         0.8, 0.2)
for (i in 1: length(allVertexName)){
  pattern = paste("^", allVertexName[i], "$", sep="")
  if (length(grep(pattern,hubRegs100)) > 0 ){
    V(g)$label[i] <- allVertexName[i]
  }else{
    V(g)$label[i] = ""
  }
}

E(g)$color <- rep("lightgray",length(E(g)))

par(mar=c(0,0,0,0))
plot(g,vertex.size = V(g)$size,
     edge.color = E(g)$color,
     vertex.label=V(g)$label,
     vertex.label.cex  = V(g)$label.cex,
     vertex.label.color = "black",
     vertex.color = V(g)$color,
     vertex.frame.color =  V(g)$color,
     #      layout = layout.reingold.tilford
     layout=layout.fruchterman.reingold(g, niter=5000, area=30*vcount(g)^2)
)

dev.off()
