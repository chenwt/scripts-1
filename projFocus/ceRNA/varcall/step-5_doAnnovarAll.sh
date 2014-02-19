#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -l h_vmem=8G,time=2::
uname -a
msg="<MESSAGE>: "
cmd="<COMMAND>: "
warn="<WARNING>: "
err="<ERROR>: "	
echo "$msg Start $0: `date`"
	
DIR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/"

##This script converts a  vcf file with one or more samples into a format required by annovar and then fires annovar's tool to get the variants.
#INPUT - $1 - a VCF file containing
#OUTPUT - $1.annovar (the input require for annovar)
#	- two other files with the filtered and dropped variants.
echo "Begin Annovar Preprocessing : Input VCF file $1";
$DIR/convert2annovar.pl $1  -format vcf4 -includeinfo -allallele > $1.annovar
echo "Converted VCF to Annovar compatible file : $1.annovar";

echo "Begin Annovar Annotation";
$DIR/summarize_annovar1.pl  $1.annovar $DIR"/humandb/" -outfile $1.annovar.summary -ver1000g 1000g2012feb -verdbsnp 135 -genetype=refgene --buildver hg19

exitstatus=$?
if [[ $exitstatus == 0 ]]
then

	echo "Finished Annovar Annotation : $1.annovar.summary.*";

	# Convert annovar result back to VCF using original VCF for header
	echo "Converting Annovar Summary to VCF : $1.summary.exome/genome_summary.csv";
	perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/convert_annovar_vcf-all-samples.pl $1.annovar.summary.exome_summary.csv $1 /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/VCFheader_forAnnovar.txt 
	perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/convert_annovar_vcf-all-samples.pl $1.annovar.summary.genome_summary.csv $1 /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/VCFheader_forAnnovar.txt 
	echo "Annovar COMPLETE ";
	echo "OUTPUT - $1.annovar.summary.exome_summary.csv.vcf ";
	echo "OUTPUT - $1.annovar.summary.genome_summary.csv.vcf ";

	rm -rf $1.annovar.summary.hg19_*
	rm -rf $1.annovar.summary.exonic_variant_function
	rm -rf $1.annovar.summary.variant_function
	rm -rf $1.annovar
	rm -rf $1.annovar.summary.invalid_input
	rm -rf $1.annovar.summary.log 
#	/ifs/scratch/c2b2/ngs_lab/sz2317/bin/src/tabix-0.2.5/bgzip -c $1.annovar.summary.genome_summary.csv.vcf > $1.annovar.summary.genome_summary.csv.vcf.gz
#	/ifs/scratch/c2b2/ngs_lab/sz2317/bin/src/tabix-0.2.5/tabix $1.annovar.summary.genome_summary.csv.vcf.gz
#	/ifs/scratch/c2b2/ngs_lab/sz2317/bin/src/tabix-0.2.5/bgzip -c $1.annovar.summary.exome_summary.csv.vcf > $1.annovar.summary.exome_summary.csv.vcf.gz
#	/ifs/scratch/c2b2/ngs_lab/sz2317/bin/src/tabix-0.2.5/tabix $1.annovar.summary.exome_summary.csv.vcf.gz

	echo "$msg End $0: `date`"
	# exit 0
fi

echo "$msg End $0: `date`"
exit $exitstatus

