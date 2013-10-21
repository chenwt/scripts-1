#!/bin/bash
#$ -cwd
#input:<full path of vcf file which was annotated by gatk and annovar old version> <-E/-W>
#output: <filtered data>
i=$1
ExomeOrWGS="-W"
GLOBAL="/ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/global_setting_GRC37lite.sh"

. $GLOBAL
	##Get Germline Variant stats - filter out variants based on MQ0Frac and Qual

	# #Filter Somatic Mutations
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/filter_vcf_somatic.rb  -v $i --dp4 2 --dp4indel 5 --dp 10 --pl 40 --plindel 60 --clr 20 --normal 2 --goesp 0.01 > $i.summary 
	if [[ $? != 0 ]]
	then
		echo "ERROR: Step 1 Filtering $i "
		exit 1
	fi
	echo "Step 1 Filtering Somatic variants .. Done"
	# ##Annotate further the filtered file with COSMIC  
	/nfs/apps/python/2.6.5/bin/python /ifs/home/c2b2/ac_lab/jh3283/tools/mutable/mutableannotateVCF.py -c /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/cancergenecensus.txt -v /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicCodingMuts_v65_28052013_noLimit.vcf -d /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/drugNameIDGene.csv -m /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicMutantExportIncFus_v66_250713.tsv -i $i.filtered.vcf> $i.filtered.vcf.ann.vcf
	if [[ $? != 0 ]]
	then
		echo "ERROR: Step 2 Cosmic Annotation $i "
		exit 1
	fi
	echo "Step 2 Cosmic Annotation .. Done"
	#convert VCF to TSV
	# ##Add extra column indicating if the variant came from exome / wgs sample's calling
	if [[ $ExomeOrWGS == "-E" ]]
	then
		extracolumn="Exome"
	else
		extracolumn="WGS"
	fi
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/convert_vcf_tsv.rb $1.filtered.vcf.ann.vcf 2 $extracolumn > $i.filtered.vcf.ann.vcf.tsv
	if [[ $? != 0 ]]
	then
		echo "ERROR:Step 3 Converting VCF to TSV $i "
		exit 1
	fi	
	echo "Step 3 Converting VCF to VCF .. Done"
	# #Add more stats by looking at the tsv
	echo "Coding Somatic Het calls that PASS filtering" >> $i.summary
	grep PASS $i.filtered.vcf.ann.vcf.tsv | grep HETEROZYGOUS |cut -f9 | sort  | uniq -c | sed 's/\s\+//' >> $i.summary

	##Annotate further the tsv with Drugbank + Clinical Trials
	/nfs/apps/python/2.6.5/bin/python /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/mutableAddClinicalTrialsToTSV.py -c  /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/ClinicalTrials.csv -d /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/DrugNameIDSynonymsTargetgeneTargetSynonyms.csv -i $i.filtered.vcf.ann.vcf.tsv > $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv
	if [[ $? != 0 ]]
	then
		echo "ERROR: Step 4 Clinical Trials and Drugbank Annotation $i "
		exit 1
	fi	
	echo "Step 4 Clinical Trials and Drugbank Annotation .. Done"
	if [[ -f FINAL.$i ]]
	then 
		rm $i.FINAL
	fi
	ln -s $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv $i.FINAL
# 	##Prioritization : 
# 	CNVFILE=`ls $Project/$WGS/$CNVDir/$wgsSampleName"-"*"-"final/$wgsSampleName"somatic_corrected_NoModLong-3-type.seg"`
# 	echo "CNVFILE is $CNVFILE"
# 	if [[ ! -e $CNVFILE || ! -s $CNVFILE || $CNVFILE == "" ]]
# 	then
# 		echo "ERROR: Step 5 Prioritization with CNVs  - $CNVFILE cnv intervals file doesnot exist"
# 	else
# 		/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/R/R-2.15.0/bin/R CMD BATCH --slave "--args $CNVFILE $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv.prioritised.tsv " /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/MutablePrioritize.R $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv.prioritised.log 
# 		if [[ $? != 0 ]]
# 		then
# 			echo "ERROR: Step 5 Prioritization with CNVs $i "
# 			doexit 1
# 		fi	
# 		echo "Step 5 Prioritization with CNVs .. Done"
# 	fi

# echo "Filtering with strict filters .. Done"