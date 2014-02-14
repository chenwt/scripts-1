#!/bin/bash
#By: J.He
#TODO: 

##in CNV reprot folder
wd=/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/CNV/report
cd $wd

# grep PARGVC TN.xcnv.cnv.bed.geneAnnot.pidtgene |sort|uniq > temp-PARGVC-TN.txt
# grep PARGVC RN.xcnv.cnv.bed.geneAnnot.pidtgene |sort|uniq > temp-PARGVC-RN.txt

# grep -v -f temp-PARGVC-TN.txt temp-PARGVC-RN.txt > temp-PARGV-RN-TN.txt
##filter CNVs at first
##use Q_Start Q_END cutoff 20

filterXcnv() {
  # cut=10
  infile=$1
  outfile=$2
  awk 'NR==1{print $0}' $infile > $outfile
  awk 'BEGIN{FS="\t"}NR>1&&$12>10&&$13>10{print $0}' $infile >> $outfile
}

# filterXcnv TN.xcnv TN.xcnv.filtered
# filterXcnv RN.xcnv RN.xcnv.filtered
grep PARGVC *N.xcnv.filtered
echo "##-----END-------"

