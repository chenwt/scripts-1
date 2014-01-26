#!/usr/bin/Rscript
#J.HE
#input: 
#output: 
#Descripiton: exploring the clinical information of current selected patients 

# Path: /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/plotSampleInfo.r


sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/sampleInfo")
}else if(sysInfo['sysname']=="Linux" ){
  print("loaded from Linux")
  setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/sampleInfo")
}

cwd=getwd()
print(cwd)

data = read.table("TCGA_barcode_all_in_cnv_meth_snp_EXP_clinical.csv",header=T)

