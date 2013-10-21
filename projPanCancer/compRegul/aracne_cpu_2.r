 ##--------------------------------------
 # get arguments from command line
 args <- commandArgs(TRUE)
 if (is.null(args)){
   print("Please provide parameters")
   exit
 }else{
   print(args)
 }
 ##--------------------------------------
tumAcro <- args[1]
num <- args[2]
type <- args[3]
pandir <- Sys.getenv("PANDIR")
tumAcro <- "blca"
num <-"2" # specify either first half or second half
type <- "sig" # specify either tf or sig

setwd(paste(pandir,tumAcro,sep=""))

expset <- paste(tumAcro,"_tcga_rnaseq_random_half",num,sep="")
output <- paste(expset,type,sep="")

##--------------------------------------
# Variables ####################
if (type == "tf") {
  regulators <- "/ifs/data/c2b2/ac_lab/malvarez/databases/tfs/tfcotf-homo-current"
} else if (type == "sig"){
  regulators <- "/ifs/data/c2b2/ac_lab/malvarez/databases/signalome/signalome-homo-current"
} else {
  print ("wrong type!(tf or sig")
  exit
}
pval<-1e-8
dpi<-0
no_bootstraps <- 100
gene<-FALSE
correct_name<-TRUE

##--------------------------------------
#print the parameters
print(expset)
print(output)
print( pval)
print( no_bootstraps)
print( regulators)
##--------------------------------------

# cheking input
if (expset == "") {
  expset <- list.files(pattern=".exp")
  if (length(expset) != 1) stop("No expset is present", call.=F)
  expset <- output <- sub(".exp", "", expset[1])
}
library(marina)



# Config threshold
if (paste("config_threshold_", output, ".txt", sep="") %in% dir()) file.rename(paste("config_threshold_", output, ".txt", sep=""), "config_threshold.txt")
cc <- 0
while (!("config_threshold.txt" %in% dir()) & cc<3) {
  cat("addpath(\'/ifs/scratch/c2b2/ac_lab/malvarez/aracne/scripts\');\nd=importdata(\'", paste(expset, ".exp", sep=""), "\');\ngenerate_mutual_threshold_configuration(d.data, \'adaptive_partitioning\');\n", sep="", file="mip.m")
  system("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfs/apps/matlab/current/runtime/glnxa64:/nfs/apps/matlab/current/bin/glnxa64:/nfs/apps/matlab/current/sys/os/glnxa64:/nfs/apps/matlab/current/sys/java/jre/glnxa64/jre/lib/amd64/server")
  system("export MCR_CACHE_ROOT=$TMPDIR/$count")
  system("export MATLABROOT=/nfs/apps/matlab/current")
  system("/nfs/apps/matlab/current/bin/matlab -singleCompThread -nodisplay -nodesktop -nojvm -nosplash < mip.m >& mip.log")
  cc <- cc + 1
}
if (!("config_threshold.txt" %in% dir())) {
  stop("Failed after 3 attempts to compute the config_threshold file", call.=F)
}

# Regulators
if (!(paste(expset, ".tfs", sep="") %in% dir())) {
  tfs <- sapply(strsplit(readLines(paste(regulators, ".txt", sep="")), "\t"), function(x) x[1])
  annot <- t(sapply(strsplit(readLines(paste(expset, ".exp", sep=""))[-1], "\t"), function(x) x[1:2]))
  probes <- annot[sub("g", "", annot[, 2]) %in% tfs, 1]
  probes <- c(probes, annot[which(!(annot[, 1] %in% probes))[1], 1])
  cat(probes, sep="\n", file=paste(expset, ".tfs", sep=""))
}
rm(tfs, annot, probes)

# ARACNe
if (!(paste(output, "-4col.txt", sep="") %in% dir())) {
  if (!(output %in% dir())) dir.create(output)
  fn <- list.files(path=paste("./", output, sep=""), pattern=".adj")
  if (length(fn)<no_bootstraps) {
    if (length(fn)>0) file.rename(paste("./", output, "/", fn, sep=""), paste("./", output, "/a", fn, sep=""))
    cat("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfs/apps/matlab/current/runtime/glnxa64
	    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfs/apps/matlab/current/bin/glnxa64
	    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfs/apps/matlab/current/sys/os/glnxa64
	    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfs/apps/matlab/current/sys/java/jre/glnxa64/jre/lib/amd64/server
	    export MCR_CACHE_ROOT=/ifs/scratch/c2b2/ac_lab/fmg2117/aracnetmp/$count
	    export MATLABROOT=/nfs/apps/matlab/current
	    /nfs/apps/matlab/current/bin/matlab -singleCompThread -nodisplay -nodesktop -nojvm -nosplash < running_aracne.m >& matlab1.log\n", sep="", file="submit_aracne.sh")
    cat("tname = \'", substr(expset, 1, 4), "\';
	    src_dir = \'/ifs/scratch/c2b2/ac_lab/malvarez/aracne/scripts\';
	    addpath(src_dir)
	    src_dir = addTrailingSlash(src_dir);
	    aracne_path = \'/ifs/scratch/c2b2/ac_lab/malvarez/aracne/ARACNE2/ARACNE\';
	    filename_exp = \'", expset, ".exp\';
	    filename_tf = \'", expset, ".tfs\';
	    method = 'adaptive_partitioning';
	    pvalue = ", pval, "; 
	    dpi = ", dpi, ";
	    no_bootstraps = ", no_bootstraps-length(fn), ";
	    bootstrap_dir = \'", output, "\';
	    verbose = 1;
	    template = 'submit-template-titan.txt';
	    template_consensusnetwork = 'consensus_template.pl';
	    aracne(template,template_consensusnetwork,src_dir,aracne_path,filename_exp,filename_tf,method,pvalue,dpi,no_bootstraps,bootstrap_dir,verbose,tname)\n", sep="", file="running_aracne.m")
    system(paste("sh ", getwd(), "/submit_aracne.sh", sep=""))
    cat("\n\nPlease run this script again when all aracne jobs were executed.\n")
    stop("", call.=F)
  }
  # Consensus and regulon
  aracne.getConsensus(output, paste(expset, ".exp", sep=""), paste(expset, ".tfs", sep=""))
}


regul <- aracne2regulon(paste(output, "-4col.txt", sep=""), paste(expset, ".exp", sep=""), "config_threshold.txt", gene=gene, mi.method="titan")
if (correct_name) {names(regul) <- sub("g", "", names(regul)); regul <- lapply(regul, function(x) {names(x$tfmode) <- sub("g", "", names(x$tfmode)); x}); names(regul) <- sub("p", "", names(regul)); regul <- lapply(regul, function(x) {names(x$tfmode) <- sub("p", "", names(x$tfmode)); x})}
class(regul) <- "regulon"
save(regul, file=paste(output, "-regulon.rda", sep=""))

# Cleaning-up
file.rename("config_threshold.txt", paste("config_threshold_", output, ".txt", sep=""))
file.remove("mip.m", "getMIThreshold.m", "mip.log", "aracnegpu.sh")
system(paste("rm -r ", output, sep=""))
