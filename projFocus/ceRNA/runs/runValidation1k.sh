#!/bin/bash
#by: j.he
#todo: 

###------functions--------
extractsnp(){
  file=$1
  awk '/^rs/&&$2>0{print $1}' $file >> temp_snp.$file
  echo "#done--$file"
}

###-------end---functions-----




#----run
###link coeff result to current wd
# while read line
# do
#   echo $line
#   ln -s ../../result/grpreg/test/grplasso_${line}.txt .
# done < /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test/gcGenes_GeneVarNet_0.05.genelist.txt 
  
# for f in `ls grplasso_*txt`
# do 
#   extractsnp $f
# done
# wc -l temp_snp.* >> GCgene_snp.counts.txt
# echo "#DONE------"


# awk '$2!="total"&&NR>4{gsub("temp_snp.grplasso_coeff","",$2);gsub(".txt","",$2);;print $2}' GCgene_snp.counts.txt > gene.list
# while read g
# do
#   extractsnp grplasso_coeff$g.txt
#   grep -wf temp_snp.grplasso_coeff$g.txt ALL.2of4intersection.20100804.sites.rsID.txt > temp_snp.grplasso_coeff$g.rsid.1kg 
# done < gene.list

# python=~/tools/python/Python_current/python
# while read gene
# do
#   $python ~/scripts/projFocus/ceRNA/validSNP.py -i grplasso_coeff${gene}.txt -d ALL.2of4intersection.20100804.sites.rsID.txt -o ${gene}_rsid_new.txt
# done < gene.list

extractsnp grplasso_coeffRANBP9.txt 
grep -wf temp_snp.grplasso_coeffRANBP9.txt ALL.2of4intersection.20100804.sites.rsID.txt > temp_snp.grplasso_coeff$g.rsid.1kg 
