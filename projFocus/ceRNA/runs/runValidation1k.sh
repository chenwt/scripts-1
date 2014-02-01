#!/bin/bash
#by: j.he
#todo: 

###------functions--------
extractsnp(){
  file=$1
  awk '/^rs/&&$2>0{print $0}' $file >> temp_snp.$file
  echo "#done--$file"
}

###-------end---functions-----


#----run
# ln -s ../../result/grpreg/test/gcgenes_genevarnet_0.05.genelist.txt gcgenes_genevarnet_0.05.genelist.txt
# ln -s ../../result/grpreg/test/grplasso_coeff_grplasso_ccnd1.txt res_ccnd1.txt
# ln -s ../../result/grpreg/test/grplasso_coeff_grplasso_gdi2.txt res_gdi2.txt
# ln -s ../../result/grpreg/test/grplasso_coeff_grplasso_nt5c1b.txt res_nt5c1b.txt
for f in `ls res_*txt`
do 
  extractsnp $f
done
wc -l temp_snp.res* >> GCgene_snp.counts.txt
rm temp*
echo "#DONE------"
