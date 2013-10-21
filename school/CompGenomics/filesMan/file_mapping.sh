

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

cat LAML_WGS.txt | awk '{print "WGS/"$1"/"$24,$16}'| sort > CHECK/files_want.txt
find WGS/ -type f -ls |grep bam |awk '{print $11}' | sort > CHECK/files_downloaded.txt
join -j1 -o 2.1 1.2 <(sort files_want.txt) <(sort files_downloaded.txt) | sort -k 2 > files_mapped.txt

# file and supposed size
cat LAML_WGS.txt | awk '{print "$23,WGS/"$1"/"$24}' | tail -n +2 