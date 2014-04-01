##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep.26,13 
#Last Updated: Mar. 1, 2014
#Usage: plot Preppi network
#Description:given gene list and presnet preppi network
# help(igraph.plotting)
source("~/scripts/myR/general.r")

args = getArgs()
cwd = system("pwd",intern=T)
require(igraph)
file =  args[1]

filenet = paste(file,".preppi",sep="")
data1 = read.table(filenet)
data2 = graph.data.frame(data1,directed=FALSE)
g = simplify(data2)
orderedVertex = get.data.frame(g,what = "vertices")
g = simplify(data2)
seedGene = unlist(read.table(file, stringsAsFactors=F))

V(g)$color <- ifelse(vapply(orderedVertex$name,FUN=function(x){length(grep(x,seedGene))},1) > 0,
                     "red", "lightblue")
V(g)$size <- ifelse(vapply(orderedVertex$name,FUN=function(x){length(grep(x,seedGene))},1) > 0,
                     5, 2)
V(g)$label.cex <- ifelse(vapply(orderedVertex$name,FUN=function(x){length(grep(x,seedGene))},1) > 0,
                    0.5, 0.2)
V(g)$label =  NA
V(g)$label = vapply(orderedVertex$name,FUN=function(x){seedGene[match(x,seedGene)]},"a") 

pdf(paste(cwd,"/plotNetwork_",tail(unlist(strsplit(filenet,"/")),1),"_label.pdf",sep=""))
par(mar=c(0,0,0,0))
plot(g,vertex.size = V(g)$size,
     vertex.label=V(g)$label,
     vertex.label.cex  = V(g)$label.cex,
     vertex.label.color = "black",
     vertex.color = V(g)$color,
     vertex.frame.color =  V(g)$color,
     layout=layout.fruchterman.reingold)
dev.off()
