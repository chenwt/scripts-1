#!/bin/bash 
#J.HE
#Input: <tab delimted file, oldfilename <tab> newfilename>
#Ouput: rename all listed files using new filename 
#PURPOSE: 

count=0
while read LINE
do	
	new=$(echo ${LINE} | awk '{print $2}')
	old=$(echo ${LINE} | awk '{print $1}')
	if [ ! -f $new ]
	then
	  cmd="mv $old $new"
	  $cmd
	  let count++
	else
	  echo "file $new exist!"
	fi	  

done < $1
echo "$count files reanmed"
echo -e "#----------------END------------------"
