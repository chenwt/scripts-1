#!/bin/bash
#$ -cwd
#By: J.He
#Desp.: prepare expression(tumor, normal), cnv, somatic(promoter, utr3, utr5), methylation matrix(diff)
#TODO: 

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/
CDT=`date|awk '{print $2"-"$3"-"$6}'`

####------get expression matrix
# resultExpDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp

# expTumDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/rnaseq/tumor
# expTumMatrix=$resultExpDir/brca_exp_l3_tumor_${CDT}.matrix 
# getExpMatfromRnaseqLevel3 $expTumDir $expTumMatrix 

# expNormDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/rnaseq/normal
# expNormMatrix=$resultExpDir/brca_exp_l3_normal_${CDT}.matrix
# getExpMatfromRnaseqLevel3 $expNormDir $expNormMatrix 

## run normalization
# Rscript step1-2_expVoomNormed.r --tumor <tumor.mat file> --normal <normal.mat file> 

# ####----get sample names
# resultGslistDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist
# expTumMatrix=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix
# head -1 $expTumMatrix | tr "\t" "\n" | awk 'NR>1'  > $resultGslistDir/barcode_tumor_Mar-21-2014.list

# ####----- get CNV matrix
# cnvTumDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/cnv/tumor
# resultCnvDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv
# genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/CG_allGene.txt
# cnvTumMatrix=$resultCnvDir/brca_cnv_l3_tumor_${CDT}.matrix
# # getCNVMat $cnvTumDir $genelist $cnvTumMatrix &
# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/getGeneSingleCnv.py -i brca_cnv_l3_tumor_Mar-23-2014.matrix -o brca_cnv_l3_tumor_Mar-23-2014.matrix.uniq.matrix


###------get methy matrix
methTumDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/methy/tumor
resultMethDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/CG_allGene.txt
methTumMatrix=$resultMethDir/brca_meth_l3_tumor_${CDT}.matrix
###1051 patients
# genMethMat $methTumDir $genelist $methTumMatrix  & #### this version is toooooooooo slow 
# genMethMat_v2 $methTumDir $genelist $methTumMatrix  &

methNormDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/methy/normal
resultMethDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/CG_allGene.txt
methNormMatrix=$resultMethDir/brca_meth_l3_normal_${CDT}.matrix
####123 patients
# genMethMat $methNormDir $genelist $methNormMatrix  &
genMethMat_v2 $methNormDir $genelist $methNormMatrix  &
##get 27number
# grep -n `grep HumanMethylation27 /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/methy/normal/input_linkfile.txt.temp|tail -1 |awk -F/ '{print $NF}'` /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/methy/normal/input_linkfile.txt.temp 
## 27
# grep -n `grep HumanMethylation27 /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/methy/tumor/input_linkfile.txt.temp|tail -1 |awk -F/ '{print $NF}'` /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/methy/tumor/input_linkfile.txt.temp
### 315

##get methydiff.matrix
# run Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/step1-5_getDiffMethy.r --tumor $methTumMatrix --normal $methNormMatrix --normal --numT 315 --numN 27

echo "#-----END---"
