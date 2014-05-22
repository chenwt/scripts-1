##after running regression to indentify candidate driver ceRNA, plot heatmap of target-regulator
## to show the predictibility of ceRNA regulators
###----under development----
rm(list=ls())
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


##--package and functions---
source(paste(rootd, "/scripts/myR/general.r",sep=""))



# args = getArgs()
# cwd = system("pwd",intern=T)
# require(igraph)
# file =  args[1]

netfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_flex_max_1000.tsv.significant.summary.netfile"

mutTFbsfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/mut_in_TFBindSite_05142014.hg19.uniq"
mutMirbsfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/mut_in_MiRNABindSite_05152014.hg19.uniq"

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

##--load data----
mutTF = loadMutGene(mutTFbsfile)
mutMir = loadMutGene(mutMirbsfile)

##----
intersect(mutTF$mutTF,mutMir$mutTF)

netRaw = loadNetFromKeyregsumfile(netfile)
netRaw = netRaw[-1,]

allDrRegs = unique(netRaw[,2])
length(allDrRegs)
c(drRegInTF= length(intersect(allDrRegs,mutTF$mutTF)), drRegInMir=length(intersect(allDrRegs,mutMir$mutTF)))
intersect(intersect(allDrRegs,mutTF$mutTF), intersect(allDrRegs,mutMir$mutTF))

write.table(file=paste(netfile,"_full.2col_",cdt,sep=""), netRaw, sep="\t", quote=F,row.names=F)
##-----

###----

data2 = graph.data.frame(netRaw,directed=FALSE)
g = simplify(data2)
orderedVertex = get.data.frame(g,what = "vertices")
g = simplify(data2)

targets = unique(netRaw[,1])
mutTFgene = unique(vapply(mutTF[,1], as.character,'a'))
mutMirgene = unique(vapply(mutMir[,1], as.character, 'a'))
mutAll  = unique(vapply(mutTF[,1],as.character,'a'), vapply(mutMir[,1], as.character, 'a'))

##----reg mutated network----
net2 = unique(rbind(netRaw[netRaw[,1] %in% mutAll, ], netRaw[netRaw[,2] %in% mutAll,]))
write.table(file=paste(netfile,"_BS_mut.2col_",cdt,sep=""), net2, sep="\t", quote=F,row.names=F)

dim(unique(netRaw)); dim(unique(net2))
c(regMutatedCeRNAdriver=length(unique(net2[,2])), target=length(unique(net2[,1]))) 
net2 <- net2[order(net2[,1]),]
tgene_regMutNum <- as.matrix(sort(table(net2[,1]),decreasing=T))
colnames(tgene_regMutNum) <- "BSmutRegNum"
write.table(file=paste(netfile,"_tgene_regMutNum_",cdt,sep=""), tgene_regMutNum,sep="\t",quote=FALSE)

net2 <- net2[order(net2[,2]),]
regMut_tgeneNum <- as.matrix(sort(table(net2[,2]),decreasing=T))
colnames(regMut_tgeneNum) <- "targetNum"
write.table(file=paste(netfile,"_regMut_tgeneNum_",cdt,sep=""), regMut_tgeneNum,sep="\t",quote=FALSE)
###--------------------------

##---network analysis
tarMutRegCount <- table(netRaw[netRaw[,2] %in% mutTFgene,1])
table(
      tarMutRegCount[order(tarMutRegCount,decreasing=T)]
      )

netMutRegMir <- table(netRaw[netRaw[,2] %in% mutMirgene,1])
table(netMutRegMir[order(netMutRegMir,decreasing=T)])
netRaw[netRaw[,2] %in% mutTFgene,][grep('ZEB2' , netRaw[netRaw[,2] %in% mutTFgene,1]),2]

# V(g)$color <- ifelse(vapply(orderedVertex$name,FUN=function(x){length(grep(paste("^",x,"$",sep=""),mutAll))},1) > 0,
#                      "red", "lightblue")


##--visulization--
data2 = graph.data.frame(net2,directed=FALSE)
g = simplify(data2)
orderedVertex = get.data.frame(g,what = "vertices")
g = simplify(data2)

##--network properties--
degNet2 <- degree(g, v = V(g), mode = "all")
coord3D <- layout.fruchterman.reingold(gD, dim = 3)
gbtws <- betweenness(g)

##--hub analysis---
pdf(paste(figd,"/step4-3_driverCeRNAnet_regBinding_mutnet_analysi", cdt,".pdf",sep=""))
ghub <- hub.score(g,scale=T)
plot(sort(ghub$vector), pch=16)

plot(1:length(ghub$vector),type="n", main = "hub score of sub network")
text(x=1:length(ghub$vector), y = 1000 * sort(ghub$vector), labels=names(sort(ghub$vector)), cex=0.8, font=2)

gaut <- authority.score(g, scale=T)
plot(sort(gaut$vector), pch=16)
plot(1:length(gaut$vector),type="n", main = "aut score of sub network")
text(x=1:length(gaut$vector), y = 1000 * sort(gaut$vector), labels=names(sort(gaut$vector)), cex=0.8, font=2)
dev.off()

allnode.hubscore <- as.matrix(sort(ghub$vector,decreasing=T))

intersect(rownames(allnode.hubscore)[1:100], rownames(tgene_regMutNum)[1:100])
intersect(rownames(allnode.hubscore)[1:100], rownames(regMut_tgeneNum)[1:100])

write.table(file=paste(netfile,".mutSubnet.hubscore_",cdt,sep=""), as.matrix(sort(ghub$vector,decreasing=T)),sep="\t",quote=FALSE)

##----top hub regulator sub network

##----hub ceRNA dirver network----
###---ceRNA dirver Summary----
hist(table(netRaw[,2][-1]))
hubRegs100 <- row.names(as.matrix(sort(table(netRaw[,2][-1]),decreasing=T)[1:100]))
write.table(file=paste(netfile,"_reg_tgeneNum_",cdt,sep=""), as.matrix(sort(table(netRaw[,2][-1]),decreasing=T)),sep="\t",quote=FALSE)

as.matrix(sort(table(netRaw[,1][-1]),decreasing=T)[1:100])


barplot(sort(table(netRaw[,2][-1]),decreasing=T)[1:100])

net3 = netRaw[netRaw[,2] %in% hubRegs100, ]
length(unique(net3[,1])); length(unique(net3[,2]))

data2 <- graph.data.frame(as.data.frame(net3),directed=FALSE)
g <- simplify(data2)

orderedVertex = get.data.frame(g,what = "vertices")
allVertexName <- orderedVertex$name

str(V(g))
for (i in 1: length(allVertexName)){
  pattern = paste("^", allVertexName[i], "$", sep="")
  if (length(grep(pattern,mutTFgene)) > 0 & length(grep(pattern, mutMirgene)) > 0 ){
    V(g)$color[i] <- "red"
    V(g)$label[i] <- allVertexName[i]
  }else if(length(grep(pattern,mutTFgene)) == 0 & length(grep(pattern, mutMirgene)) > 0){
    
    V(g)$color[i] <- "orange"
    V(g)$label[i] <- allVertexName[i]
    
  }else if(length(grep(pattern,mutTFgene)) > 0 & length(grep(pattern,mutMirgene)) == 0) {
    
    V(g)$color[i] <- "steelblue"
    V(g)$label[i] <- allVertexName[i]
    
  }else{
    V(g)$color[i] = "gray"
    V(g)$label[i] = ""
    
  }
}

V(g)$color <- ifelse(vapply(orderedVertex$name,FUN=function(x){pattern = paste("^", x, "$", sep=""); length(grep(pattern,hubRegs100))},1) > 0,
                    "red", "lightblue")
V(g)$size <- ifelse(vapply(orderedVertex$name,FUN=function(x){pattern = paste("^", x, "$", sep="");length(grep(pattern,hubRegs100))},1) > 0,
                    2, 0.8)

V(g)$label.cex <- ifelse(vapply(orderedVertex$name,FUN=function(x){pattern = paste("^", x, "$", sep="");length(grep(pattern,hubRegs100))},1) > 0,
                         0.5, 0.1)
for (i in 1: length(allVertexName)){
  pattern = paste("^", allVertexName[i], "$", sep="")
  if (length(grep(pattern,hubRegs100)) > 0 ){
    V(g)$label[i] <- allVertexName[i]
   }else{
    V(g)$label[i] = ""
  }
}

# pdf(paste(figd,"/step4-3_hubAnalsyis_hubReg_mutnet_",cdt,".pdf",sep=""))
pdf(paste(figd,"/step4-3_hubAnalsyis_hubReg_",cdt,".pdf",sep=""),width=10,)
E(g)$color <- rep("lightgray",length(E(g)))
par(mar=c(0,0,0,0))
plot(g,vertex.size = V(g)$size,
     edge.color = E(g)$color,
     vertex.label=V(g)$label,
     vertex.label.cex  = V(g)$label.cex,
     vertex.label.color = "black",
     vertex.color = V(g)$color,
     vertex.frame.color =  V(g)$color,
     layout = layout.drl 
#      layout = layout.reingold.tilford
#      layout=layout.fruchterman.reingold(g, niter=3000, area=30*vcount(g)^2)
)
dev.off()

par(mar=c(5,5,1,1))
max(degree(g))
pdf(paste(figd,"/step3-4_degreeDist_test.pdf",sep=""))
plot(degree.distribution(g,v=V(g)),xlim=c(1,117),xlab="node degree")

lines(degree.distribution(g))
dev.off()

##---community detection---
lec <- leading.eigenvector.community(g)

##---visulization--
orderedVertex = get.data.frame(g,what = "vertices")
allVertexName <- orderedVertex$name
V(g)$label =  NA
for (i in 1: length(allVertexName)){
  pattern = paste("^", allVertexName[i], "$", sep="")
  if (length(grep(pattern,mutTFgene)) > 0 & length(grep(pattern, mutMirgene)) > 0 ){
    V(g)$color [i] <- "red"
    V(g)$label[i] <- allVertexName[i]
  }else if(length(grep(pattern,mutTFgene)) == 0 & length(grep(pattern, mutMirgene)) > 0){
    V(g)$color [i] <- "orange"
    V(g)$label[i] <- allVertexName[i]
    
  }else if(length(grep(pattern,mutTFgene)) > 0 & length(grep(pattern,mutMirgene)) == 0) {
    V(g)$color [i] <- "steelblue"
    V(g)$label[i] <- allVertexName[i]
    
  }else{
    V(g)$color [i] <- "gray"
  }
}

V(g)$size <- ifelse(vapply(orderedVertex$name,FUN=function(x){length(grep(x,targets))},1) > 0,
                    2, 0.8)

# V(g)$label =  NA
# V(g)$label = vapply(orderedVertex$name,FUN=function(x){targets[match(x,targets)]},"a") 

pdf(paste(figd,"/figure_4-3_optCoorMutatedNetork)",cdt,".pdf",sep=""), width = 10, height= 10)
par(mar=c(0,0,0,0))
plot(g,vertex.size = V(g)$size,
     vertex.label=V(g)$label,
     vertex.label.cex  = V(g)$label.cex,
     vertex.label.color = "black",
     vertex.color = V(g)$color,
     vertex.frame.color =  V(g)$color,
     layout=layout.fruchterman.reingold(g, niter=3000, area=30*vcount(g)^2)
    )
dev.off()
