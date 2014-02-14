#!/bin/bash
#By: J.He
#TODO: 



samtools mpileup -ugf ref.fa aln.bam | bcftools view -bvcg - > var.raw.bcf 
bcftools view var.raw.bcf | vcfutils.pl varFilter -D 100 > var.flt.vcf

