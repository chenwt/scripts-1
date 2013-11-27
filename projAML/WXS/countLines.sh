#!/bin/bash
#By: J.He
#TODO: 
#input:<file: list of full path of file names to count lines>
#output: <file: list of lines in the listed files>

outfile=$1".countlines"
for file in `cat $1`
do
  wc -l $file >> $outfile
done   
echo "#----------END----"

