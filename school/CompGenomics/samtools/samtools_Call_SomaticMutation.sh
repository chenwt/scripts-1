
 ###------------------------------call raw variants
MYREF='/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/samtools/wu_build36.fasta'
# if [ ! -f "samtools_Call_SomaticMutation.sh" ] ;then  cp ../samtools_Call_SomaticMutation.sh . ; fi
#samtools mpileup -ugSDf $MYREF tumor.$CHROM.bam normal.$CHROM.bam | bcftools view -bvcgT pair -p 1.1 - > $CHROM.var.bcf
#echo "Pileup $CHROM Done!" >> logsrun
bcftools view $CHROM.var.bcf | vcfutils.pl varFilter -D 300 > $CHROM.var.vcf

#call somatic variantsq
ruby ~yshen/scratch/do_filter-somatic.rb -v $CHROM.var.vcf
echo "Calling somatic variants $CHROM DONE!" >> logsrun

