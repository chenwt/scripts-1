#!/bin/bash
#$ -cwd
#input:<full path of vcf file which was annotated by gatk and annovar old version> <-E/-W>
#output: <filtered data>

GLOBAL=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh
. $GLOBAL
i=$1
	##Get Germline Variant stats - filter out variants based on MQ0Frac and Qual
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/filter_vcf_coding.rb -v $i
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/vcf_transition-transversion-per-sample.rb $i.germline-filtered.vcf 1 > $i.summary.germline.coding.stats
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/vcf_transition-transversion-per-sample.rb $i.germline-filtered.vcf -1 > $i.summary.germline.noncoding.stats
	# rm -rf $i.germline-filtered.vcf $i.germline-dropped.vcf

	#Filter Somatic Mutations
	# $RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/filter_vcf_somatic.rb  -v $i --dp4 2 --dp4indel 5 --dp 10 --pl 40 --plindel 60 --clr 20 --normal 2 --goesp 0.01 > $i.summary 
	echo "filter somatic..."
	$RUBY /ifs/home/c2b2/ac_lab/jh3283/tools/mutable/filter_vcf_somatic_Jing.rb -v $i --dp4 2 --dp4indel 5 --dp 10 --pl 40 --plindel 60 --clr 20 --normal 2 --goesp 0.01 > $i.summary 

	if [[ $? != 0 ]]
	then
	  echo "ERROR: ruby filtering somatic Region Step 1 Filtering $i "
		exit 1
	fi
	# echo "Step 1 Filtering Somatic variants .. Done"
	##Annotate further the filtered file with COSMIC  
	## New python script version
	# python /ifs/home/c2b2/ac_lab/jh3283/tools/mutable/mutableannotateVCF_Jing.py -c /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/cancergenecensus.txt -v /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicCodingMuts_v65_28052013_noLimit.vcf -d /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/drugNameIDGene.csv -m /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicMutantExportIncFus_v66_250713.tsv -i $i.germline-filtered.vcf > $i.germline-filtered.vcf.ann.vcf
	# if [[ $? != 0 ]]
	# then
	# 	echo "ERROR: python Step 2 Cosmic Annotation $i "
	# 	exit 1
	# fi
	# echo " Step 2 Cosmic Annotation .. Done"

	# #convert VCF to TSV
	# ##Add extra column indicating if the variant came from exome / wgs sample's calling
      	# extracolumn="WGS"
	# $RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/convert_vcf_tsv.rb $1.germline-filtered.vcf.ann.vcf 2 $extracolumn > $i.germline-filtered.vcf.ann.vcf.tsv
	# if [[ $? != 0 ]]
	# then
	# 	echo "ERROR: Exonic Region Step 3 Converting VCF to TSV $i "
	# 	exit 1
	# fi	
	# echo "Exonic Region Step 3 Converting VCF to TSV .. Done"

  	
   	echo "Filtering  Cancer Mode with strict filters .. Done"
