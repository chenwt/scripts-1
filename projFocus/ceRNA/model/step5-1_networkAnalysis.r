# 
# install.packages("sna")
# install.packages("statnet");update.packages("statnet")
# install.packages("rgl")
# install.packages("coda")
# install.packages("numDeriv")
# install.packages("yacca")
# install.packages("informR")
# demo(package="igraph")

##after running regression to indentify candidate driver ceRNA, plot heatmap of target-regulator
## to show the predictibility of ceRNA regulators
###----under development----
rm(list=ls())
require("igraph");require(sna)
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

# netfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_flex_max_1000.tsv.significant.summary.netfile"
# netfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/keyRegSummary_allRunMay_05212014_0.01"
netfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_flex_max_1000.tsv.significant.summary.netfile.cancertargetgene"

mutTFbsfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/mut_in_TFBindSite_05142014.hg19.uniq"
mutMirbsfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/mut_in_MiRNABindSite_05152014.hg19.uniq"
expfile = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"

###---func---
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
netRaw = loadNetFromKeyregsumfile(netfile)
netRaw = netRaw[-1,]

data2 = graph.data.frame(netRaw,directed=FALSE)
g = simplify(data2)
orderedVertex = get.data.frame(g,what = "vertices")

expD <- read.table(expfile, header=T,sep="\t")
exp2D <- na.omit(expD[unlist(orderedVertex), ])
summary(exp2D)
rm(expD)

targets = unique(netRaw[,1])
mutTFgene = unique(vapply(mutTF[,1], as.character,'a'))
mutMirgene = unique(vapply(mutMir[,1], as.character, 'a'))
mutAll  = unique(c(vapply(mutTF[,1],as.character,'a'), vapply(mutMir[,1], as.character, 'a')))

require(tsne)
exp2D_tsne <- tsne(exp2D,k=3,max_iter=100)
require(rgl)
plot3d(exp2D_tsne[rownames(exp2D) %in% mutAll == FALSE, ],col = "gray", xlab="x", ylab="y",zlab="z" )
plot3d(exp2D_tsne[rownames(exp2D) %in% mutAll == TRUE,], col = "red",add = T)

plot(exp2D_tsne[,1:2])
points(exp2D_tsne[rownames(exp2D) %in% mutAll,1:2], col = "red")
# exp2D_pca <- princomp(exp2D)


##----reg mutated network----

data2 = graph.data.frame(net2,directed=FALSE)
g = simplify(data2)
orderedVertex = get.data.frame(g,what = "vertices")

data2 = graph.data.frame(netRaw,directed=FALSE)
g = simplify(data2)
orderedVertex = get.data.frame(g,what = "vertices")

g = simplify(data2)

allVertexName <- orderedVertex$name
for (i in 1: length(allVertexName)){
  pattern = paste("^", allVertexName[i], "$", sep="")
  if (length(grep(pattern,mutTFgene)) > 0 & length(grep(pattern, mutMirgene)) > 0 ){
    V(g)$color [i] <- "red"
    V(g)$label[i] <- allVertexName[i]
  }else if(length(grep(pattern,mutTFgene)) == 0 & length(grep(pattern, mutMirgene)) > 0){
    V(g)$color [i] <- "orange"
    V(g)$label[i] <- allVertexName[i]
    
  }else if(length(grep(pattern,mutTFgene)) > 0 & length(grep(pattern,mutMirgene)) == 0) {
    V(g)$color [i] <- "green"
    V(g)$label[i] <- allVertexName[i]
    
  }else{
    V(g)$color[i] <- "lightblue"
    V(g)$label[i] <- ''
  }
}

V(g)$color <- ifelse(vapply(allVertexName,FUN=function(x){pattern = paste("^", x, "$", sep="");length(grep(pattern,mutAll))},1) > 0, 'red', 'lightblue')

V(g)$size <- ifelse(vapply(allVertexName,FUN=function(x){pattern = paste("^", x, "$", sep="");length(grep(pattern,mutAll))},1) > 0, 3, 2)

V(g)$label.cex <- ifelse(vapply(allVertexName,FUN=function(x){pattern = paste("^", x, "$", sep="");length(grep(pattern,mutAll))},1) > 0, 0.6, 0.1)

E(g)$color <- rep("lightgray",length(E(g)))
pdf(paste(figd,"/sep4-3_summary_cancerGene_optCorr_mutatedNetowrk.", cdt, ".2.pdf",sep=""))
par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(g,vertex.size = V(g)$size,
     edge.color = E(g)$color,
     vertex.label=V(g)$label,
     vertex.label.cex  = V(g)$label.cex,
     vertex.label.color = "black",
     vertex.color = V(g)$color,
     vertex.frame.color =  V(g)$color,
     #      layout = layout.drl 
     #      layout = layout.reingold.tilford
     layout=layout.fruchterman.reingold(g, niter=3000, area=30*vcount(g)^2)
)
dev.off()

# tkplot(g, layout=layout.kamada.kawai)
l <- layout.reingold.tilford

l <- layout.drl(g)
rglplot(g, layout=l)

summary(g)

##---analysis---
cumty_ws <- walktrap.community(g)
cumty_fg <- fastgreedy.community(g)
cumty_lev <- leading.eigenvector.community(g)
cumty_im <- infomap.community(g)
par(mfrow=c(2,2))
plot(cumty_ws, g)
plot(cumty_fg, g)
plot(cumty_lev, g)
plot(cumty_im, g)

modularity(cumty_fg); modularity(cumty_ws); modularity(cumty_im); modularity(cumty_lev)
barplot(table(membership(cumty_fg)[mutAll]))
barplot(table(membership(cumty_fg)))

kruskal.test(table(membership(cumty_fg)), table(membership(cumty_fg)[mutAll]))



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

V(g)$size <- ifelse(vapply(orderedVertex$name,FUN=function(x){pattern = paste("^", x, "$", sep="");length(grep(pattern,hubRegs100))},1) > 0, 5, 1)

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
par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(g,vertex.size = V(g)$size,
     edge.color = E(g)$color,
     vertex.label=V(g)$label,
     vertex.label.cex  = V(g)$label.cex,
     vertex.label.color = "black",
     vertex.color = V(g)$color,
     vertex.frame.color =  V(g)$color,
#      layout = layout.drl 
     #      layout = layout.reingold.tilford
          layout=layout.fruchterman.reingold(g, niter=3000, area=30*vcount(g)^2)
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
