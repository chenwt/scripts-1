#!/bin/bash
for BAMLIST in `ls *.list`
do
BAM=`echo $BAMLIST |sed 's/\@//g'` 
echo "#!/bin/bash" >> $BAM.sh
echo "BAM="$BAM >> $BAM.sh
#echo GATKJAR=/ifs/data/c2b2/ngs_lab/ngs/usr/GATK/Sting/dist/GenomeAnalysisTK.jar
#echo GATK="java -Xmx3072m -jar ""$GATKJAR"
#echo regionlist=/ifs/scratch/c2b2/ac_lab/jh3283/ref/Exome_Targeted_Regions.BED
#echo REF=/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta
#echo outputname=$BAM

echo java -Xmx3072m -jar /ifs/data/c2b2/ngs_lab/ngs/usr/GATK/Sting/dist/GenomeAnalysisTK.jar \
 -T DepthOfCoverage \
 -I $BAM \
 -L /ifs/scratch/c2b2/ac_lab/jh3283/ref/Exome_Targeted_Regions.BED \
 -R /ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta \
 -dt BY_SAMPLE \
 -dcov 5000 \
 -l INFO \
 -o $BAM \
 --omitDepthOutputAtEachBase \
 --omitLocusTable \
 --minBaseQuality 0 \
 --minMappingQuality 20 \
 --start 1 \
 --stop 5000 \
 --nBins 200 >> $BAM.sh
qsub -l mem=4G,time=2:: -N $BAM -cwd $BAM.sh
done
