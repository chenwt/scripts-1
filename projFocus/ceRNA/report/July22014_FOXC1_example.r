rm(list=ls())

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
source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

figd = paste(rootd,"/DATA/projFocus/report/Jul2014/fig",sep="")
fcnet=paste(rootd,"/DATA/projFocus/report/Jul2014/net.cernet",sep="")
flnet=paste(rootd,"/DATA/projFocus/report/Jul2014/net.lassonet",sep="")
flnetmut=paste(rootd, "/DATA/projFocus/report/Jul2014/net.lassonet.mutated",sep="")
fgnet=paste(rootd,"/DATA/projFocus/report/Jul2014/net.greedynet",sep="")


require(igraph)

plotNet = function(fcnet, verCol) {
  cnet = read.table(fcnet,sep="\t")
  net = graph.data.frame(cnet)
  g = simplify(net)
  orderedVertex = get.data.frame(g,what = "vertices")
  ###plot out properties
  coreGene = 'FOXC1'
  vcolor <- vapply(orderedVertex$name, FUN=function(x){
                  ifelse(x==coreGene,'red',verCol)},'a')
  vsize <- vapply(orderedVertex$name, FUN=function(x){
                  ifelse(x==coreGene,8,6)},1.0)
  vlabel = vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,x,"")},'a')
  
  vlabelsize <- vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,1,0.5)},1.0)
  ecolor <- rep("lightgray",length(E(g)))
  
  ##--plot
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  plot(g,
       vertex.size = vsize,
       edge.color = ecolor,
       vertex.label = vlabel,
       vertex.label.cex  = vlabelsize,
       vertex.label.color = "black",
       vertex.color = vcolor,
       vertex.frame.color =  vcolor,
       edge.arrow.size=0.1,
       layout=layout.fruchterman.reingold(g, niter=300, area=30*vcount(g)^2)
  )

}

plotCandiNet = function(fgnet, mutgene) {
  gnet = read.table(fgnet,sep="\t")
  net = graph.data.frame(gnet)
  g = simplify(net)
  orderedVertex = get.data.frame(g,what = "vertices")
  ###plot out properties
  coreGene = 'FOXC1'
  vcolor <- vapply(orderedVertex$name, FUN=function(x){
    pattern = paste("^",x,"$",sep="")
    if(x==coreGene) {
      return("red")
    }else if(length(grep(pattern, mutgene)) ==1) {
      return("lightblue")
    }else {
      return("gray")
    }},
    'a')
  
  vsize <- vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,8,6)},1.0)
  vlabel = vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,x,'')},'a')
  
  vlabelsize <- vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,1.3,0.7)},1.0)
  ecolor <- rep("lightgray",length(E(g)))
  
  ##--plot
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  plot(g,
       vertex.size = vsize,
       edge.color = ecolor,
       vertex.label = vlabel,
       vertex.label.cex  = vlabelsize,
       vertex.label.color = "black",
       vertex.color = vcolor,
       vertex.frame.color =  vcolor,
       edge.arrow.size=0.1,
       layout=layout.fruchterman.reingold(g, niter=300, area=30*vcount(g)^2)
  )
  
}

plotMutNet = function(fgnet, mutgene) {
  gnet = read.table(fgnet,sep="\t")
  net = graph.data.frame(gnet)
  g = simplify(net)
  orderedVertex = get.data.frame(g,what = "vertices")
  ###plot out properties
  coreGene = 'FOXC1'
  vcolor <- vapply(orderedVertex$name, FUN=function(x){
    pattern = paste("^",x,"$",sep="")
    if(x==coreGene) {
      return("red")
    }else if(length(grep(pattern, mutgene)) ==1) {
      return(colors()[411])
    }else {
      return("lightblue")
    }},
    'a')
  
  vsize <- vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,8,6)},1.0)
  vlabel = vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,x,'')},'a')
  
  vlabelsize <- vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,1.3,0.7)},1.0)
  ecolor <- rep("lightgray",length(E(g)))
  
  ##--plot
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  plot(g,
       vertex.size = vsize,
       edge.color = ecolor,
       vertex.label = vlabel,
       vertex.label.cex  = vlabelsize,
       vertex.label.color = "black",
       vertex.color = vcolor,
       vertex.frame.color =  vcolor,
       edge.arrow.size=0.1,
       layout=layout.fruchterman.reingold(g, niter=300, area=30*vcount(g)^2)
  )
  
}

plotMutNet2 = function(flnet, mutgene, selectgene) {
  lnet = read.table(flnet,sep="\t")
  net = graph.data.frame(lnet)
  g = simplify(net)
  orderedVertex = get.data.frame(g,what = "vertices")
  ###plot out properties
  coreGene = 'FOXC1'
  mutgene = setdiff(mutgene, selectgene)
  vcolor <- vapply(orderedVertex$name, FUN=function(x){
    pattern = paste("^",x,"$",sep="")
    if(x==coreGene) {
      return("red")
    }else if(length(grep(pattern, mutgene)) ==1) {
      return(colors()[411])
    }else if(length(grep(pattern, selectgene)) ==1) {
      return(colors()[499])
    }else {
      return("lightblue")
    }},
    'a')
  

  vsize <- vapply(orderedVertex$name, FUN=function(x){
    ifelse(x==coreGene,8,6)},1.0)
  
  vlabel <- vapply(orderedVertex$name, FUN=function(x){
    pattern = paste("^",x,"$",sep="")
    if(x==coreGene) {
      return(x)
    }else if(length(grep(pattern, mutgene)) ==1) {
      return("")
    }else if(length(grep(pattern, selectgene))==1){
      return(x)    
    }else {
      return("")
    }},
    'a')
  
  vlabelsize <- vapply(orderedVertex$name, FUN=function(x){
    pattern = paste("^",x,"$",sep="")
    if(x==coreGene) {
      return(1.3)
    }else if(length(grep(pattern, mutgene)) ==1) {
      return(0.5)
    }else if(length(grep(pattern, selectgene))==1){
      return(0.9)    
    }else {
      return(0.5)
    }},
    1.2)
  ecolor <- rep("lightgray",length(E(g)))
  
  ##--plot
  par(mar=c(0,0,0,0),mfrow=c(1,1))
  plot(g,
       vertex.size = vsize,
       edge.color = ecolor,
       vertex.label = vlabel,
       vertex.shape = "circle",
       vertex.label.cex  = vlabelsize,
       vertex.label.color = "black",
       vertex.color = vcolor,
       vertex.frame.color =  vcolor,
       edge.arrow.size=0.1,
       layout=layout.fruchterman.reingold(g, niter=300, area=30*vcount(g)^2)
  )
  
}


lasgene = read.table(flnet,sep="\t",stringsAsFactors=F)[,2]
mutgene = read.table(flnetmut,stringsAsFactors=F)[,1]
gredgene=read.table(fgnet,sep="\t",stringsAsFactors=F)[,2]

pdf(paste(figd,"/report_July_example_FOXC1.pdf",sep=""))
plotNet(fcnet,'gray')
plotCandiNet(fcnet,lasgene)
plotNet(flnet, "lightblue") 
plotMutNet(flnet, mutgene)

plotMutNet2(flnet, mutgene, gredgene)
dev.off()
