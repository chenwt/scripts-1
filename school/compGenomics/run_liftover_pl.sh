## run liftover.pl

## source for chain file http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/liftOver/


VCF='/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_132.b37.vcf'
CHAIN='/ifs/scratch/c2b2/ys_lab/jh3283/ref/hg19ToHg18.over.chain'
OUT='/ifs/scratch/c2b2/ys_lab/jh3283/ref/dbsnp_130.b36.vcf'
GATK='/ifs/data/c2b2/ngs_lab/ngs/resources/GATK-2.3-9/'
NEWREF='/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/wu_build36.fasta'
OLDREF='/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/Homo_sapiens_assembly19.fasta'

 perl /ifs/home/c2b2/ys_lab/jh3283/scripts/school/compGenomics/lifeOverVCF.pl -vcf $VCF -chain $CHAIN -out $OUT -gatk $GATK -newRef $NEWREF -oldRef $OLDREF