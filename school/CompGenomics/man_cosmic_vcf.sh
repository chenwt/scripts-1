# manipulate cosmic vcf file

egrep -v "^#|^$" hg19_cosmic_v54_furgason.vcf | sed 's/chr//' >  hg19_cosmic_v54_furgason.vcf.new

grep "^#" hg19_cosmic_v54_furgason.vcf > header="header.txt"

cat hg19_cosmic_v54_furgason_new.vcf >> header.txt
mv header.txt hg19_cosmic_v54_furgason_new.vcf