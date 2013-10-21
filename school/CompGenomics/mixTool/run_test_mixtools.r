
if(Sys.info()["sysname"] == "Darwin") {
	setwd("/Volumes/ys_lab_scratch/jh3283/school/compGenomic/mixtools/wd/")
	print("work on IMAC!")
}else {
	setwd("~/SCRATCH/school/compGenomic/mixtools/wd/")
	print("Work On TITAN!")
}


source("/ifs/home/c2b2/ys_lab/jh3283/scripts/school/CompGenomics/mixTool/test_mixtools.r")

# fnames <- dir()[grep("freq",dir())]
fnames <- "2966.somaticfiltered.freq"
data <- getData(fnames)


##------------------------------CV to select lambda for mixture model
#one sample
kmix <- 10
baf <- data$BAF
loglikes <- myMultiMixMLE(baf,kmix)
kMixture <- getKMix(baf,loglikes)
plotNormalComponent(kMixture,length(kMixture$lambda))
#multiple samples
#input: files ,output: lambda, mu, list
nf <- length(fnames)
for (i in 1: nf){
	

}


test
lambda <- 0.5; mu <- c(0.3, 0.7) 
mix1 <- myGaussianEM(getData(fnames)$BAF,lambda,mu)

out <- "res_2966"
# sapply(seq(.1,.9,by=.1),function(x){myGaussianEM(data$BAF,x,mu)$lambda})

loglike.normalmix(data$BAF,mixture=mix1)




	