
rm(list=ls())
rdafile = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_tcga_rnaseq851_signalome-regulon.rda"

load(rdafile)
output = paste(rdafile,".4col.txt",sep="")

outheader = c("tf","tfmod","val1","likelihood")
write.table(file=output, paste(outheader,collapse="\t"),col.names=F,quote=F,row.names=F,sep="\t")

for (i in 1:length(regul)){
  tflikhd_crt = as.matrix(regul[[i]]$likelihood)
  ntfmod = length(tflikhd_crt)
  out_crtD = as.data.frame(matrix(NA, ncol=4, nrow=ntfmod))
  tfmode_crt = regul[[i]]$tfmode
  out_crtD[,1] = rep(names(regul[i]), ntfmod)
  out_crtD[,2] = names(tfmode_crt)
  out_crtD[,3] = tfmode_crt
  out_crtD[,4] = regul[[i]]$likelihood
  write.table(file=output,out_crtD,col.names=F,row.names=F,quote=F,sep="\t",append=T)
}


rm(list=ls())
rdafile = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_tcga_rnaseq851_geneid-regulon.rda"
load(rdafile)
output = paste(rdafile,".4col.txt",sep="")

outheader = c("tf","tfmod","val1","likelihood")
write.table(file=output, paste(outheader,collapse="\t"),col.names=F,quote=F,row.names=F,sep="\t")
for (i in 1:length(regul)){
  tflikhd_crt = as.matrix(regul[[i]]$likelihood)
  ntfmod = length(tflikhd_crt)
  out_crtD = as.data.frame(matrix(NA, ncol=4, nrow=ntfmod))
  tfmode_crt = regul[[i]]$tfmode
  out_crtD[,1] = rep(names(regul[i]), ntfmod)
  out_crtD[,2] = names(tfmode_crt)
  out_crtD[,3] = tfmode_crt
  out_crtD[,4] = regul[[i]]$likelihood
  write.table(file=output,out_crtD,col.names=F,row.names=F,quote=F,sep="\t",append=T)
}
