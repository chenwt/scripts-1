#!/bin/bash
#Author: Jing He
#Date: Apr.17
#Last Update: 
#Example:  
DISABR="LAML"
LBS="exome"

tail -n +2 LAML_exome.txt | awk '{print "exome/"$1"/"$24,$23,$16}'| sort > CHECK/files_LAML_exome_want.txt
find exome/ -type f -ls |grep bam |awk '{print $11}' | sort > CHECK/files_LAML_exome_downloaded.txt
join -j1 -o 2.1 1.2 <(sort CHECK/files_LAML_exome_want.txt) <(sort CHECK/files_LAML_exome_downloaded.txt) | sort -k 2 > CHECK/files_LAML_exome_mapped.txt

grep -f CHECK/files_LAML_exome_mapped.txt CHECK/files_LAML_exome_want.txt > CHECK/files_LAML_sucess_exome.txt
cat files_LAML_sucess_exome.txt | sed -e 's/\//\t/g' | awk '{if($5 ==11) print $1"/"$2"/"$3 }' > files_LAML_sucess_exome_Normal.txt
cat files_LAML_sucess_exome.txt | sed -e 's/\//\t/g' | awk '{if($5 ==03) print $1"/"$2"/"$3 }' > files_LAML_sucess_exome_Tumor.txt


cat files_LAML_sucess_exome_Normal.txt | sed -e 's/\//\t/g' | cut -f3 | cut -d- -f3 > normalId.sucess.exome.temp
grep -f normalId.sucess.exome.temp files_LAML_sucess_exome_Tumor.txt  | sed 's/\//\t/g' | sort -k 3 | awk '{print $1"/"$2"/"$3}' > files_exome_Paired_tumor.txt

cut -d/ -f3 files_exome_Paired_tumor.txt | cut -d- -f3 > tumorId.sucess.exome.temp
grep -f tumorId.sucess.exome.temp files_LAML_sucess_exome_Normal.txt | sed 's/\//\t/g' | sort -k 3 | awk '{print $1"/"$2"/"$3}' > files_exome_Paired_normal.txt


