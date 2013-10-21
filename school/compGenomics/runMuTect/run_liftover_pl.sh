## run liftover.pl

## source for chain file http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/liftOver/
## example: ~/scripts/school/compGenomics/runMuTect/run_liftover_pl.sh

VCF='/ifs/scratch/c2b2/ys_lab/jh3283/ref/test/hg19_cosmic_v54_furgason.vcf'
CHAIN='/ifs/scratch/c2b2/ys_lab/jh3283/ref/hg19ToHg18.over.chain'
OUT='/ifs/scratch/c2b2/ys_lab/jh3283/ref/test/hg18_cosmic_v54_furgason.vcf'
GATK='/ifs/home/c2b2/ys_lab/jh3283/tools/GenomeAnalysisTK-2.4-9-g532efad'
# NEWREF='/ifs/scratch/c2b2/ys_lab/jh3283/ref/Homo_sapiens_assembly18.fasta'
OLDREF='/ifs/scratch/c2b2/ys_lab/jh3283/ref/Homo_sapiens_assembly18'
# OLDREF='/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/wu_build36'
# OLDREF='/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/Homo_sapiens_assembly19'
NEWREF='/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/Homo_sapiens_assembly19'

# echo $VCF
# echo $CHAIN
# echo $GATK
# echo $OUT


cmd="perl  /ifs/home/c2b2/ys_lab/jh3283/scripts/school/compGenomics/runMuTect/lifeOverVCF.pl -vcf $VCF -chain $CHAIN -out $OUT -gatk $GATK -newRef $NEWREF -oldRef $OLDREF"
echo $cmd
$cmd