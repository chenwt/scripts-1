###set global config
export REFTYPE="build"  ## build: no "chr" in chromosome names, aka >1 >2 etc;  hg: >chr1 >chr2 etc.
export DBSNP="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_130_b37.rod"
export DBSNPVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_132.b37.excluding_sites_after_129.vcf" 
export DBSNP132="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_132.b37.vcf"
export DBSNP135="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_135_b37.vcf"
export HapMapV3VCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/hapmap_3.3.b37.vcf"
export INDELVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/1000G_indels_for_realignment.b37.vcf"
export OneKGenomes="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/ALL.wgs.phase1.projectConsensus.snps.sites.vcf.PASS.vcf"
export EVSVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/EVS.vcf"
##we need this for gatk VQSR new version
export OMNI_Sites="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GATK_resourceBundle/1000G_omni2.5.b37.sites.vcf"
export HamMap_Sites="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GATK_resourceBundle/hapmap_3.3.b37.sites.vcf"
export Mills1KG_Sites="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GATK_resourceBundle/Mills_and_1000G_gold_standard.indels.b37.sites.vcf"

export SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
export BCFTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/bcftools"
export GATKJAR_old="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK/Sting/dist/GenomeAnalysisTK.jar"
export GATKJAR12="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK-1.2-1/GenomeAnalysisTK.jar"
export GATKJAR="/ifs/home/c2b2/ac_lab/jh3283/tools/GATK/GenomeAnalysisTK_current/GenomeAnalysisTK.jar"
export ANNOVAR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/"	##this is exome veriosn. chk if match wid wgs version
export PICARD="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools-1.96/"   # export PICARD="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools-1.65/"
export STING="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK"
export AnnotationTable="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/refGene-big-table-b37.txt"
export RUBY18="/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/ruby"
export RUBY="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/ruby"
export UTILS="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/"
export NGSSHELL="/ifs/data/c2b2/ngs_lab/ngs/code/shell/"
export FIXMATE="$PICARD/FixMateInformation.jar"
export FIXSTAT="$PICARD/CollectInsertSizeMetrics.jar"
export GCbias="$PICARD/CollectGcBiasMetrics.jar"


## obsolete since WGS is also now using BWA mem , hence the new bwa-0.7.3 version and new indexes of the ref genome
## BWA-aln settings for WGS
## export REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
## export BWA="/ifs/data/c2b2/ngs_lab/ngs/Pipeline/bwa-0.5.9"
## export BPATH="/ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/"

## BWA-mem settings for Exome & WGS
export REF="/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa"
export BWA="/ifs/data/c2b2/ngs_lab/ngs/usr/src/bwa-0.7.3/bwa"
export BPATH="/ifs/home/c2b2/ac_lab/jh3283/tools/mutable"

# For exome seq the default targets are Illumina Nextera.
export ExonFile="/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/TruSeq-Exome-Targeted-Regions-BED-file.sorted.bed"
export GenesTarget="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/AllGenes-b37.bed"


export ListRuns="/RunsPerSample/"
export SAMPLES="SAMPLES/"
export INFO="INFO/"
export somaticCalling="somaticCalling/"
export germlineCalling="germlineCalling/"
export CNVDir="CNV/hmmCopy/"
export LOG="/LOG/"
export WGS="/WGS/"
export Exome="/Exome/"
export RNA="/RNA/"
export RESULTS="/Results/"


