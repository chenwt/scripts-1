#!/bin/bash 
#J.HE
#Input: <tab delimted file, oldfilename <tab> newfilename>
#Ouput: rename all listed files using new filename 
#PURPOSE: 

count=0
log=${0##*/}.log
echo -n "" >$log
while read LINE
do	
	new=$(echo ${LINE} | awk '{print $2}')
	old=$(echo ${LINE} | awk '{print $1}')
	if [ ! -f $new ]
	then
	  cmd="ln -s $old $new"
	  $cmd
	  echo $new >>$log
	  let count++
	fi	  

done < $1
echo "$count files reanmed"
echo -e "#----------------END------------------"
