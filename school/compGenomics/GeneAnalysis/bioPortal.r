#!/bin/Rscirpts
#Author: Jing He
#Date: May 7
#Last Updat: 
#parameters: 
#example: Rscript 

###------------------------------header------------------------------
args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters")
	exit
}else{
	print(args)
}

###------------------------------start coding##------------------------------

if(Sys.info()["sysname"] == "Darwin"){
		Sys.setenv(COMPGENOMDIR="/Volumes/ys_lab_scratch/jh3283/school/compGenomic"
				   ,SCRIPTS="/Users/jh3283/Dropbox/scripts/"
	      )
	}
COMPGENOMDIR <- Sys.getenv("COMPGENOMDIR")
SCRIPTS <- Sys.getenv("SCRIPTS")
PIDFILE <- paste(COMPGENOMDIR,"/samtools/dirList.txt",sep="")
setwd(COMPGENOMDIR)
###------------------------------start coding##------------------------------

install.packages("cgdsr")
require(cgdsr)
myPID <- unlist(read.table(PIDFILE,header=F))

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[1,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[4,1]

# Get data slices for a specified list of genes, genetic profile and case list
getProfileData(mycgds,c('BRCA1','BRCA2'),mygeneticprofile,mycaselist)


# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist[3,1])

mycaseClinicalData <- myclinicaldata[gsub("TCGA.AB.","",rownames(myclinicaldata)) %in% 
c("2966",as.character(as.vector(myPID))),]

write.table(mycaseClinicalData,"patients_Clinical_Infor.txt",quote=FALSE,sep="\t")

# documentation
help('cgdsr')
help('CGDS')