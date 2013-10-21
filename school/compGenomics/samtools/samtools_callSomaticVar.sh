samtools mpileup -DSuf $MYREF tumor.$1.bam normal.$1.bam | bcftools view -bvcgT pair -p 1.1 - > $1.var.bcf
bcftools view $1.var.bcf | vcfutils.pl varFilter -D 60 > $1.var.vcf
