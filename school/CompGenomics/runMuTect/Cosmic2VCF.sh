# #!/bin/sh -cwd 
#Author: Jing He
#Date: Mar. 29, 2013
#Last Update: 
#Example:  Cosmic2VCF.sh file

cat $1 | cut -f12,18 > pre.$1 
tail -n +2 pre.$1 | awk '{print $1}' | sed 's/.*\(...\)$/\1/'| sed 's/>/\t/g' > f1.temp

tail -n +2 $1 | awk '{print $2}' | sed 's/\:\|\-/\t/g'| awk '{print $1"\t"$3}' > f2.ChrPos.temp
## file format CHROM POS ID<.> REF ALT QUAL<.> FILTER<.> INFO<SOMATIC>
tail -n +2 $1 | awk '{print $1}' | sed 's/.*/\./g'  >f3.dot.temp
tail -n +2 $1 | awk '{print $1}' | sed 's/.*/SOMATIC/g'  >f4.INFO.temp
paste f2.ChrPos.temp f3.dot.temp f1.bb.temp f3.dot.temp f3.dot.temp f4.INFO.temp  > $1.temp
cat $1.temp | awk '{ if ($4 ~ /^[ATGC]/ && $5 ~ /^[ATGC]/) print $0 }' > $1.temp.temp

cat hg19_cosmic_v54_furgason.vcf | grep "^#" > VCFHead.file.temp
cat VCFHead.file.temp $1.temp > $1.out
# #head sed 's/hg19/b36/g'

# rm *temp
