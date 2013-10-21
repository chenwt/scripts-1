if [ -f $1.BED ] ;then  rm $1.BED ; fi
if [ -f $1.vcf ] ;then  rm $1.vcf ; fi

cd $1
cat 1.var.vcf.somaticfiltered.vcf | grep "^#" > $1.vcf
for chrom in `seq 1 22` X Y  
do
	cat $chrom.var.vcf.somaticfiltered.vcf | awk '{if (!/^#/){print $0}}' >> ../wd/$1.vcf
done

echo "$1.vcf created!"
	# sed -e 's/chr//' $1.vcf | awk '{OFS="\t"; if (!/^#/){print "chr"$1,$2-1,$2,$3,$4}}' >> $1.BED

# ~/tools/liftOver/liftOver $1.BED ~/SCRATCH/ref/hg18ToHg19.over.chain $1.hg19.BED $1.left.BED




