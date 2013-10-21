#!/bin/bash
#Author: Jing He
#Date: 
#Last Update: 
#Example:  



##------------------------------initializing
$DIS_ABBR=$1
$CWD=

#generate .txt file from xml file
Rscripts ~/scripts/school/compGenomic/parseXML.r 

cat ../../cgquery/LAML_WGS.txt | awk '{print $26,"  ",$25}'  > CHECK/wgs.md5

# cut LAML_WGS.txt -f1,13,17,26,25 | head

# 
cat LAML_WGS.txt | awk '{print "WGS/"$1"/"$24,$16}'| sort

# list all downloaded files
find WGS/ -type f -ls |grep bam |awk '{print $11}' | sort

tail -n +2 LAML_WGS.txt | awk '{print "WGS/"$1"/"$24,$16}'| sort > CHECK/files_want.txt
find WGS/ -type f -ls |grep bam |awk '{print $11}' | sort > CHECK/files_downloaded.txt
join -j1 -o 2.1 1.2 <(sort files_want.txt) <(sort files_downloaded.txt) | sort -k 2 > files_mapped.txt

grep -f files_downloaded_WGS.txt files_want_WGS.txt > files_sucess_WGS.txt
cat files_sucess_WGS.txt | sed -e 's/\//\t/g' | awk '{if($5 ==11) print $1"/"$2"/"$3 }' > Files_Sucess_WGS_SolidTissueNormal.txt
cat files_sucess_WGS.txt | sed -e 's/\//\t/g' | awk '{if($5 ==03) print $1"/"$2"/"$3 }' > Files_Sucess_WGS_PrimaryBloodDerivedCancer.txt


cat Files_Sucess_WGS_SolidTissueNormal.txt | sed -e 's/\//\t/g' | cut -f3 | cut -d- -f3 > normalId.sucessWGS.temp
grep -f normalId.sucessWGS.temp Files_Sucess_WGS_PrimaryBloodDerivedCancer.txt  | sed 's/\//\t/g' | sort -k 3 | awk '{print $1"/"$2"/"$3}' > WGS_Paired_PBDC.txt

cut -d/ -f3 WGS_Paired_PBDC.txt | cut -d- -f3 > tumorIDsucessWGS.temp
grep -f tumorIDsucessWGS.temp Files_Sucess_WGS_SolidTissueNormal.txt  | sed 's/\//\t/g' | sort -k 3 | awk '{print $1"/"$2"/"$3}' > WGS_Paired_Normal.txt

#  output pairwise sample file names
# /ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/LAML/CHECK/WGS_Paired_Normal.txt
# /ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/LAML/CHECK/WGS_Paired_PBDC.txt
##------------------------------exome  
# cat LAML_exome.txt | awk '{print "exome/"$1"/"$24,$16}'| sort

# # list all downloaded files
# find exome/ -type f -ls |grep bam |awk '{print $11}' | sort

tail -n +2 LAML_exome.txt | awk '{print "exome/"$1"/"$24,$23,$16}'| sort > CHECK/files_want.exome.txt
find exome/ -type f -ls |grep bam |awk '{print $11}' | sort > CHECK/files_exome_downloaded.txt
join -j1 -o 2.1 1.2 <(sort files_want.exome.txt) <(sort files_exome_downloaded.txt) | sort -k 2 > files_mapped.txt

grep -f files_exome_downloaded.txt files_want.exome.txt > files_sucess_exome.txt
cat files_sucess_exome.txt | sed -e 's/\//\t/g' | awk '{if($5 ==11) print $1"/"$2"/"$3 }' > Files_Sucess_exome_Normal.txt
cat files_sucess_exome.txt | sed -e 's/\//\t/g' | awk '{if($5 ==03) print $1"/"$2"/"$3 }' > Files_Sucess_exome_Tumor.txt


cat Files_Sucess_WGS_SolidTissueNormal.txt | sed -e 's/\//\t/g' | cut -f3 | cut -d- -f3 > normalId.sucessWGS.temp
grep -f normalId.sucessWGS.temp Files_Sucess_WGS_PrimaryBloodDerivedCancer.txt  | sed 's/\//\t/g' | sort -k 3 | awk '{print $1"/"$2"/"$3}' > WGS_Paired_PBDC.txt

cut -d/ -f3 WGS_Paired_PBDC.txt | cut -d- -f3 > tumorIDsucessWGS.temp
grep -f tumorIDsucessWGS.temp Files_Sucess_WGS_SolidTissueNormal.txt  | sed 's/\//\t/g' | sort -k 3 | awk '{print $1"/"$2"/"$3}' > WGS_Paired_Normal.txt

