#!/usr/bin/env python
#J.HE
#Desp:created for projFocus to annot each somatic mutated gene by gene starting and ending position.Dec 17,2013
#     this code need to be run before connect with expression data, find the promoter region somatic mutation gene
#input:<file: .mat file generated using getSomMutMatrix.py> <file:annotation file from lab> 
#output:<file: .mat.anno which add chrom, start, end, strand informationas annotation>
#TODO:




inpAnno="/ifs/scratch/c2b2/ac_lab/rs3412/no1/net/entrez_annotation_hg19.txt"

