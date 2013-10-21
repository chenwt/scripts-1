setwd("C:\\Users\\Sam Ni\\Dropbox\\G4002\\Project")
data <- read.table("data/0221_jing/exportqueryresults.csv",sep=",",stringsAsFactors=F)
colnames(data) <- c("index","subject_id","icustay_id","ioitemid","charttime","elemid","cgid","cuid","site","rate","rateuom")

require(ggplot2)
ggplot(data,aes(x=subject_id,y=icustay_id)) + geom_point(color="#ffcccc") + opts(title="NICU_subid_stayid")
 
