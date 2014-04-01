##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep.26,13 
#Last Updated: Mar. 1, 2014
#Usage: plot Preppi network
#Description:given gene list and presnet preppi network
# help(igraph.plotting)
source("~/scripts/myR/general.r")
source("~/scripts/myR/plot/colorRampPaletteAlpha/colorRampPaletteAlpha.R")
# source("~/HOME//scripts/myR/general.r")
# source("~/HOME/scripts/myR/plot/colorRampPaletteAlpha/colorRampPaletteAlpha.R")

args = getArgs()
cwd = "/ifs/data/c2b2/ac_lab/jh3283/projMisc/yLan//Mar29"
print(cwd)
require(igraph)
file1 =  args[1] 
file2 = args[2]
outfile = paste(cwd,"/plotNetwork_",tail(unlist(strsplit(file1,"/")),1),"_",
            tail(unlist(strsplit(file2,"/")),1),"_label.pdf",sep="")
print(args)
#---tests--
# library("igraph")
# setwd(("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projMisc/yLan//Mar29"))
file1="denovalGenename.list"
file2="mouseTPgenename.list"
outfile = paste('/Volumes/',cwd,"/plotNetwork_",tail(unlist(strsplit(file1,"/")),1),"_",
                tail(unlist(strsplit(file2,"/")),1),"_label.pdf",sep="")

data1 = read.table(jxy(file1,".preppi"))
data2 = read.table(jxy(file2,".preppi"))

net = graph.data.frame(rbind(data1,data2),directed=FALSE)
g = simplify(net)
orderedVertex = get.data.frame(g,what = "vertices")
g = simplify(net)
seedGene1 = unlist(read.table(file1, stringsAsFactors=F))
seedGene2 = unlist(read.table(file2, stringsAsFactors=F))
seedGene = c(seedGene1, seedGene2)
data1gene =  unlist(c(as.character(data1[,1]),as.character(data1[,2])))


layout = layout.fruchterman.reingold(g, niter=5000)
V(g)$color <- vapply(orderedVertex$name,FUN=function(x){
      if (length( grep(jxy("^",x,"$"), seedGene1) ) > 0) {return(addalpha("red",0.8))}
      else if (length(grep(jxy("^",x,"$"), seedGene2)) > 0) {return(addalpha("blue",0.5))}
      else return(addalpha("lightblue",0.7))}, 'a')
V(g)$size <- ifelse(vapply(orderedVertex$name,FUN=function(x){length(grep(jxy("^",x,"$"),seedGene))},1) > 0,
                     4,1.5) 
V(g)$label.cex <- ifelse(vapply(orderedVertex$name,FUN=function(x){length(grep(jxy("^",x,"$"),seedGene))},1) > 0,
                    0.5, 0.2)
V(g)$label =  NA
V(g)$label = sapply(orderedVertex$name,FUN=function(x){idx = grep(jxy("^", x,"$"),seedGene);
                                                       ifelse(length(idx)>0,seedGene[idx],NA) }) 

pdf(outfile)
par(mar=c(0,0,0,0))
plot(g,vertex.size = V(g)$size,
     vertex.label=V(g)$label,
     vertex.label.cex  = V(g)$label.cex,
     vertex.label.color = "black",
     vertex.label.font = 2,
     vertex.color = V(g)$color,
     vertex.frame.color =  V(g)$color,
     edge.width = 0.3,
     layout=layout)
dev.off()
