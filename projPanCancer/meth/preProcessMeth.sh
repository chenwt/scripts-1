#!/bin/bash
#By: J.He
#TODO: 

# preProcessing methylation files: # tar -zxvf; # report the number of files;  identify replicates ; getting ride of duplicates

rootd="/ifs/scratch/c2b2/ac_lab/jh3283/projPanCancer/meth/"
cd $rootd$1
d=$(ls *tar.gz |head -1)
echo $d

dtest=$(echo $d |rev |cut -b 8-|rev)
echo $dtest

#if [ ! -d ${dtest} ]; then
# for file in  $(ls *tar.gz)
# do
#   tar -zxvf ${file}
# done
#fi
#
#echo -e "Study\tTotal\tMatched\tUnmatched\n" > ${rootd}report_$1.txt
#for d in $(ls -p |grep /)
#do
#  cd ${rootd}$1"/"${d}
#  numUM=$(grep jhu MANIFEST.txt|awk '{print $2}'|cut -d- -f5|sort| uniq -c | awk '$1!=2{print $2}' | wc -l )
#  numM=$(grep jhu MANIFEST.txt|awk '{print $2}'|cut -d- -f5|sort| uniq -c | awk '$1==2{print $2}' |wc -l)
#  numTotal=$(grep jhu MANIFEST.txt|wc -l)
#  echo -e $d"\t"$numTotal"\t"$numM"\t"$numUM >> ${rootd}report_$1.txt
#  cd ${rootd}$1
#done
