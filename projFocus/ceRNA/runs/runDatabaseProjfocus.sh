#!/bin/bash
#By: J.He
#TODO: 

#cp annot_snpArray_6.0_txt.bed annot_snpArray_6.0_txt.bed.bak
#sed -i 1d annot_snpArray_6.0_txt.bed 
#sed -i "1iIdentifer\tChrom\tPos\tStrand" annot_snpArray_6.0_txt.bed
#awk '$3!="n/a"{print $0}' ucsc_gene_tss |wc 

###----------snp_annotation_file
#grep -v ^# GenomeWideSNP_6.na22.annot.csv |awk 'BEGIN{FS=",";OFS="\t";print "identifier","chromosome","pos","strand","dbsnpID"}NR>1{print $1,$4,$5,$6,$3}' |sed "s/\"//g" > annot_GenomeWideSNP_6_5cols.txt
#awk '$2!="---"{print $0}' annot_GenomeWideSNP_6_5cols.txt > annot_GenomeWideSNP_6_5cols_clean.txt
#rm annot_GenomeWideSNP_6_5cols.txt
#awk '$5!="---"{print $0}' annot_GenomeWideSNP_6_5cols_clean.txt > annot_GenomeWideSNP_6_5cols.txt
#mv annot_GenomeWideSNP_6_5cols.txt annot_GenomeWideSNP_6_5cols_clean.txt


#------------rnaseq_exp_annotation---file
#awk 'BEGIN{OFS=FS="\t";print "GeneSymbol","chrom","start","strand"}NR>1{gsub("chr","",$1);print $6,$1,$3,$2}' refseq_gene_tss_geneSymbol |uniq> refseq_gene_tss_geneSymbol.4cols 
#head -20 refseq_gene_tss_geneSymbol.4cols > test
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/get_tss_anno_byGene.py test
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/get_tss_anno_byGene.py refseq_gene_tss_geneSymbol.4cols
##---need-the-5cols-file_to_use_anno_SNP.py---
#awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,$1}' refseq_gene_tss_geneSymbol.4cols.signleTSS  > refseq_gene_tss_geneSymbol.signleTSS.5cols


#------------gene_annotation_by_start_end_pos---file
#python getSingleStartEnd.py
#python convertEntrezAnno2GeneAnno.py -i entrez_annotation_hg19_unique.txt  -o gene_annotation_hg19_unique_start_end.txt


##---------------------
#awk 'BEGIN{FS=OFS="\t"}NR>1{print "chr"$2,$3,$3+1,$4,$5":"$1}' annot_GenomeWideSNP_6_5cols_clean.txt > annot_GenomeWideSNP_NCBI36_chr.bed


##---------------
# awk 'BEGIN{FS=OFS="\t"; print "genename","chrom","strand","start","end"};
    # NR>1{gsub("chr","",$2);print $8,$2,$3,$4,$5}'  refseq_gene_tss > refset_gene_start_end.tsv 
# python ~/scripts/projFocus/ceRNA/step1-3_getGeneStartEnd.py -i refset_gene_start_end.tsv.sorted.uniq_resortedCol -o refset_gene_start_end.tsv.sorted.uniq_resortedCol_geneSingleStartEnd
