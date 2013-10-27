#!/bin/bash
#$ -cwd
#input:<full path of vcf file which was annotated by gatk and annovar old version> <-E/-W>
#output: <filtered data>
i=$1
ExomeOrWGS=$2

. /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/global_setting_GRC37lite.sh

# echo $RUBY
	##Get Germline Variant stats - filter out variants based on MQ0Frac and Qual
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/filter_vcf_coding.rb -v $i
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/vcf_transition-transversion-per-sample.rb $i.germline-filtered.vcf 1 > $i.summary.germline.coding.stats
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/vcf_transition-transversion-per-sample.rb $i.germline-filtered.vcf -1 > $i.summary.germline.noncoding.stats
	rm -rf $i.germline-filtered.vcf $i.germline-dropped.vcf

	# #Filter Somatic Mutations
	$RUBY /ifs/home/c2b2/ac_lab/jh3283/tools/mutable/filter_vcf_somatic_Jing.rb  -v $i --dp4 2 --dp4indel 5 --dp 10 --pl 40 --plindel 60 --clr 20 --normal 2 --goesp 0.01 > $i.summary 
	if [[ $? != 0 ]]
	then
		echo "ERROR: Exonic Region Step 1 Filtering $i "
		exit 1
	fi
	echo "Exonic Region Step 1 Filtering Somatic variants .. Done"
	##Annotate further the filtered file with COSMIC  
	## OLD R script /ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/R/R-2.15.0/bin/R CMD BATCH --slave "--args germline-filtered.vcf germline-filtered.vcf.ann.vcf " /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/MutableAnnotate.R germline-filtered.vcf.ann.log 
	## OLD python version
	## 	python /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/mutableannotateVCF.py -c /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/cancergenecensus.txt -v /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicCodingMuts_v65_28052013_noLimit.vcf -d /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/drugNameIDGene.csv -i germline-filtered.vcf > germline-filtered.vcf.ann.vcf
	## New python script version
	python /ifs/home/c2b2/ac_lab/jh3283/tools/mutable/mutableannotateVCF.py -c /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/cancergenecensus.txt -v /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicCodingMuts_v65_28052013_noLimit.vcf -d /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/drugNameIDGene.csv -m /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicMutantExportIncFus_v66_250713.tsv -i $i.filtered.vcf > $i.filtered.vcf.ann.vcf
	# python /ifs/home/c2b2/ac_lab/jh3283/tools/mutable/mutableannotateVCF_Jing.py -c /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/cancergenecensus.txt -v /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicCodingMuts_v65_28052013_noLimit.vcf -d /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/drugNameIDGene.csv -m /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/CosmicMutantExportIncFus_v66_250713.tsv -i $i.filtered.vcf > $i.filtered.vcf.ann.vcf
	if [[ $? != 0 ]]
	then
		echo "ERROR: Exonic Region Step 2 Cosmic Annotation $i "
		exit 1
	fi
	echo "Exonic Region Step 2 Cosmic Annotation .. Done"
	#convert VCF to TSV
	##Add extra column indicating if the variant came from exome / wgs sample's calling
	if [[ $ExomeOrWGS == "-E" ]]
	then
		extracolumn="Exome"
	else
		extracolumn="WGS"
	fi
	$RUBY /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/convert_vcf_tsv.rb $1.filtered.vcf.ann.vcf 2 $extracolumn > $i.filtered.vcf.ann.vcf.tsv
	if [[ $? != 0 ]]
	then
		echo "ERROR: Exonic Region Step 3 Converting VCF to TSV $i "
		exit 1
	fi	
	echo "Exonic Region Step 3 Converting VCF to VCF .. Done"
	#Add more stats by looking at the tsv
	echo "Coding Somatic Het calls that PASS filtering" >> $i.summary
	grep PASS $i.filtered.vcf.ann.vcf.tsv | grep HETEROZYGOUS |cut -f9 | sort  | uniq -c | sed 's/\s\+//' >> $i.summary

	##Annotate further the tsv with Drugbank + Clinical Trials
	# $PYTHON /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/mutableAddClinicalTrialsToTSV.py -c  /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/ClinicalTrials.csv -d /ifs/data/c2b2/ngs_lab/ngs/resources/mutable/DrugNameIDSynonymsTargetgeneTargetSynonyms.csv -i $i.filtered.vcf.ann.vcf.tsv > $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv
	# if [[ $? != 0 ]]
	# then
	# 	echo "ERROR: Exonic Region Step 4 Clinical Trials and Drugbank Annotation $i "
	# 	exit 1
	# fi	
	# echo "Exonic Region Step 4 Clinical Trials and Drugbank Annotation .. Done"
	# # ln -s $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv FINAL.$i

# 	##Prioritization : TODO on both Exome & WGS sample - make exome stuff wait for WGS cnv analysis
# #	CNVFILE=`ls $Project/$CNVDir/$sample"-"*"-"final/$sample"somatic_corrected_NoModLong-3-type.seg.filteredoutputfile.seg" `
# 	wgsSampleName=`cat $Project/$INFO/TUMOR `
# 	exomeSampleName=`cat $Project/$INFO/EXOMETUMOR `
# 	CNVFILE=`ls $Project/$WGS/$CNVDir/$wgsSampleName"-"*"-"final/$wgsSampleName"somatic_corrected_NoModLong-3-type.seg"`
# 	echo "CNVFILE is $CNVFILE"
# 	if [[ ! -e $CNVFILE || ! -s $CNVFILE || $CNVFILE == "" ]]
# 	then
# 		echo "ERROR: Exonic Region Step 5 Prioritization with CNVs  - $CNVFILE cnv intervals file doesnot exist"
# 	else
# 		/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/R/R-2.15.0/bin/R CMD BATCH --slave "--args $CNVFILE $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv.prioritised.tsv " /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/MutablePrioritize.R $i.filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv.prioritised.log 
# 		if [[ $? != 0 ]]
# 		then
# 			echo "ERROR: Exonic Region Step 5 Prioritization with CNVs $i "
# 			doexit 1
# 		fi	
# 		echo "Exonic Region Step 5 Prioritization with CNVs .. Done"
# 	fi


	##Copy important data to results folder & Merge the TSVs from exome sequencing and WGS sequencing analysis
	
	# ResultDir=$Project/$RESULT/
	# if [ ! -d $ResultDir ]; thefiltern mkdir $ResultDir ; fi
	# if [[ $ExomeOrWGS == "-E" ]]
	# then
	# 	cp  filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv $ResultDir/exonicvariants.ExomeAnalysis.tsv
	# 	cp  filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv.prioritised.tsv $ResultDir/exonicvariants.prioritised.ExomeAnalysis.tsv
	# 	cp $i.summary $ResultDir/exonicvariants.ExomeAnalysis.somatic.summary.stats
	# 	cp $i.summary.germline.coding.stats $ResultDir/exonicvariants.ExomeAnalysis.germline.coding.stats
	# 	cp $i.summary.germline.noncoding.stats $ResultDir/exonicvariants.ExomeAnalysis.germline.noncoding.stats			
	# else
	# 	cp  filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv $ResultDir/exonicvariants.WGSAnalysis.tsv
	# 	cp  filtered.vcf.ann.vcf.tsv.clinicaltrial.tsv.prioritised.tsv $ResultDir/exonicvariants.prioritised.WGSAnalysis.tsv
	# 	cp $i.summary $ResultDir/exonicvariants.WGSAnalysis.somatic.summary.stats
	# 	cp $i.summary.germline.coding.stats $ResultDir/exonicvariants.WGSAnalysis.germline.coding.stats
	# 	cp $i.summary.germline.noncoding.stats $ResultDir/exonicvariants.WGSAnalysis.germline.noncoding.stats			
	# fi


 #   	if [[ -e $Project/$WGS/$somaticCalling/$wgsSampleName/ALL/all.exome.vcf.filtered.vcf.ann.vcf.tsv.prioritised.tsv && -e $Project/$Exome/$somaticCalling/$exomeSampleName/ALL/all.exome.vcf.filtered.vcf.ann.vcf.tsv.prioritised.tsv ]];
 #   	then
 #   	##TODO Sort by chr-pos
 #   		cat $Project//$WGS/$somaticCalling/$wgsSampleName/ALL/all.exome.vcf.filtered.vcf.ann.vcf.tsv.prioritised.tsv <( grep -v "^Chr" $Project/$Exome/$somaticCalling/$exomeSampleName/ALL/all.exome.vcf.filtered.vcf.ann.vcf.tsv.prioritised.tsv ) >> $ResultDir/exonicvariants.ExomeAndWGSAnalysis.tsv
 #   	fi
  	
   	echo "Filtering in Personal Cancer Mode with strict filters .. Done"