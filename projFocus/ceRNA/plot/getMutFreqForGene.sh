#!/bin/bash
#! -cwd
#By: J.He
#Desp.: generate mutations/mutation group number for each gene
#TODO: 

file=$1

awk 'BEGIN{g="";n=0}{
  if (NR==1||$1==g)
    {n=n+1
    g=$1}
  else
    {print g"\t"n
    g=$1
    n=1} }
    END{print g"\t"n}
    ' $file > $file.MutFreqForGene
