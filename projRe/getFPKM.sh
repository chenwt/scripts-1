#!/bin/bash
#$ -cwd
# J.HE
total=$(ls **/genes.fpkm_tracking|wc -l)
i=1
if [ -f result_genes_fpkm_matrix_${total} ]; then rm result_genes_fpkm_matrix_${total}; fi  
for line in $(ls **/genes.fpkm_tracking)
do
    if [ ${i} != 1 ]; then
      #echo ${i}
      mv result_genes_fpkm_matrix_${total} temp
      awk 'NR>1{print $10}' ${line} | paste temp - > result_genes_fpkm_matrix_${total}
    else
      awk 'NR>1{print $1"\t"$10 }' ${line} > result_genes_fpkm_matrix_${total}
    fi
    i=$((i+1))
done

ls **/genes.fpkm_tracking | cut -d/ -f1 > barcode.list.fpkm
rm temp
echo "$((i-1)) samples"


