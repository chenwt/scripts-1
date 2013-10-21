#!/bin/bash
#$ -cwd

#Author:J.He
#Date: Jun 13,2013

# make director using tumorAcr 
# while read tumAcr ; do mkdir $tumAcr ; done < tumAcro.txt
# create softlink to exp files
#while read tumAcr 
#do
#cd ~/SCRATCH/projPanCancer/compRegul/${tumAcr}/
#rm *exp
#ln -s ~/SCRATCH/projPanCancer/PANCANCER/${tumAcr}/${tumAcr}_tcga_rnaseq_random_half1.exp ${tumAcr}_tcga_rnaseq_random_half1.exp 
#echo "linking"${tmAcr}
#ln -s ~/SCRATCH/projPanCancer/PANCANCER/${tumAcr}/${tumAcr}_tcga_rnaseq_random_half2.exp ${tumAcr}_tcga_rnaseq_random_half2.exp 
#
#done < tumAcro.txt

## run the first section of aracne code
#while read tumAcr
#do
#cd ~/SCRATCH/projPanCancer/compRegul/${tumAcr}/
#cp ../aracne_cpu_1.r .
#echo "first part of aracne run on "${tumAcro}
##Rscript aracne_cpu_1.r ${tumAcr} 1 tf 
##Rscript aracne_cpu_1.r ${tumAcr} 2 tf 
#Rscript aracne_cpu_1.r ${tumAcr} 1 sig 
#Rscript aracne_cpu_1.r ${tumAcr} 2 sig 

#done < tumAcro_input.txt
#
# moniter the aracne jobs until the last done


# run the second section of the aracne code
while read tumAcr
do
cd ~/SCRATCH/projPanCancer/compRegul/${tumAcr}/
cp ../aracne_cpu_2.r .
echo "first part of aracne run on "${tumAcro}
#Rscript aracne_cpu_2.r ${tumAcr} 1 tf 
#Rscript aracne_cpu_2.r ${tumAcr} 2 tf 
#Rscript aracne_cpu_2.r ${tumAcr} 1 sig 
Rscript aracne_cpu_2.r ${tumAcr} 2 sig 
done < tumAcro_input.txt


##--THis is the end------------------------------------

