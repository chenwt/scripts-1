#!/bin/bash
#! -cwd
#By: J.He
#Desp.: get frequency of each mutation/mutation group
#TODO: 

file=$1
echo -n "" > $file.mutFreq
while read line 
do
  # for(i=5;i<NF;i++)
  # print $1,$2,$3,$4,sum}' >> $file.mutFreq
  echo $line| awk '{
    sum=0
    for(i=8;i<NF;i++)
      {if($i!=0)
	sum=sum+1
      }
      print $1,$2,$6,$7,sum}' >> $file.mutFreq
done < $file
