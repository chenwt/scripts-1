samtools mpileup -uf wu_build36.fasta tumor.$chrom.bam | bcftools view -bvcg - > tumor.$chrom.var.raw.bcf
