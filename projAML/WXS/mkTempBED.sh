#!/bin/bash
#$ -cwd

#Author:J.He
i=1; 
while read file ;
do
  filename=$(echo $file| cut -d/ -f9 |sed "s/dbSNP135commonchr//g")
  cut -f1,2,3 $file | grep -v "^@" > temp${filename} 
done < tempBEDfiles.txt
