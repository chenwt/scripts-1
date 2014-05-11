# install.packages("plotrix")
library(plotrix)
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,4,7)],collapse="-")
figd =  "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/May2014/doc/"
output = paste(figd, "figure1_", CDT, ".pdf", sep="")
pdf(output)

parts  <- c(5352 - 3340,3340)
pct = parts/ sum(parts)

ctgs <- c("Non-dysGint", "dysGint")

mylabs <- paste(ctgs,
                paste(format(pct * 100, digits=2),rep("%",2), sep=""),
                sep = "\n")
mytitle = paste("Percentage of Gint\n(",sum(parts), ")",sep="")
mycol = colorRampPalette(c("blue","white"))(length(parts))

mycol =  c("lightblue","orange")

pie3D(parts, labels=mylabs, labelpos=lp,
      border="white",col=mycol,
      main=mytitle,
      explode=.2,radius=pi/4,shade=0.5,theta=1
      )

dev.off()