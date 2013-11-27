#!/usr/bin/Rscript
#input: <file: output from uTest_snp_exp.py> <filenames: for output>
#output: <file: with adjusted p-value,only output significant snps>

args <- commandArgs(TRUE)
 if (is.null(args)){
      print("Please provide parameters")
   exit
    }else{
         print(args)
    }
input = args[1]
outfile = paste(input,".adjusted",sep="") 
padj_cut = 10^(-6)
data_old = read.table(input ) 
colnames(data_old) = sapply(data_old[1,],as.character)
data_old = data_old[-1,] 
p_old = data_old[,4]
p_adj = p.adjust(sapply(p_old,function(x){as.numeric(as.character(x))}),method='bonferroni') 
out = cbind(data_old[which(p_adj < padj_cut),],p_adj[which(p_adj < padj_cut)])
print(nrow(out))
colnames(out) <- c(colnames(data_old),'pval_adj')
write.table(out,outfile,sep="\t",quote =F,row.names=F)
print("#--------END-------")
