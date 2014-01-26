##scripts to documentation work flow and current working status
##created for projFocus
##J.HE
##Jan.23,2014

require(R2HTML)
sysInfo = Sys.info()
ifelse(sysInfo['sysname'] == "Darwin",
       setwd("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/doc"),
       setwd("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/doc"))
cwd = getwd()
HTMLStart(outdir=cwd, file="projFocusWorkFlow",
          extension="html", echo=FALSE, HTMLframe=TRUE)
HTML.title("Project Focus ceRNA Work flow and Status", HR=1)
HTML.title("outline",HR=3)

HTMLhr()

# HTML.title("plot1", HR=2)
# plot(fresult$PrGreaterThanF)
# HTMLplot() 

HTMLStop()