###set global config
export REFTYPE="build"  ## build: no "chr" in chromosome names, aka >1 >2 etc;  hg: >chr1 >chr2 etc.
export REF="/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa"
export DBSNP="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_130_b37.rod"
export DBSNPVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_132.b37.excluding_sites_after_129.vcf" 
export DBSNP132="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_132.b37.vcf"
export HapMapV3VCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/hapmap_3.3.b37.vcf"
export INDELVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/1000G_indels_for_realignment.b37.vcf"
export OneKGenomes="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/ALL.wgs.phase1.projectConsensus.snps.sites.vcf.PASS.vcf"
export EVSVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/EVS.vcf"
export SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
export BCFTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/bcftools"
export GATKJAR="/ifs/home/c2b2/ac_lab/jh3283/tools/GATK/GenomeAnalysisTK_current/GenomeAnalysisTK.jar"
# export GATKJAR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar"
# export GATKJAR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GenomeAnalysisTK-1.5-21-g979a84a/GenomeAnalysisTK.jar"
export PICARD="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools-1.65/"
export BWA="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/bwa"  
export BPATH="/ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/WGS/Mutable/"
export STING="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK"
export AnnotationTable="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/refGene-big-table-b37.txt"
export RUBY18="/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/ruby"
export UTILS="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/"
export NGSSHELL="/ifs/data/c2b2/ngs_lab/ngs/code/shell/"
export MYSCDIR="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/"

export FIXMATE="$PICARD/FixMateInformation.jar"
export FIXSTAT="$PICARD/CollectInsertSizeMetrics.jar"
export GCbias="$PICARD/CollectGcBiasMetrics.jar"


## For exome seq the default targets are Illumina Nextera.
export ExonFile=/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/TruSeq-Exome-Targeted-Regions-BED-file.sorted.bed

