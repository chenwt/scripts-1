setwd("/Volumes/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/clinicInfo")
data = read.table("TARGET_AML_ClinicalDataSet_20120702.txt",sep="\t",header=T,fill=NA)
rownames(data) = data$USI
pid16 = readLines("PID_16.txt")
data = data[pid16,]
save(data,file="clinicInfo16Patient.rda")
summary(data)