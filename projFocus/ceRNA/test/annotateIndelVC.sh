#!/bin/bash
#$ -cwd
# J.HE
##Annotate the varaints files using GATK
GLOBAL='/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/global_setting_GRC37lite.sh'

. ${GLOBAL}
CMD="sh $BPATH/gatk_annotate_variants.sh -v $OUT -b $INP -s $GLOBAL "
echo $CMD
$CMD 1> $OUT.log.varann.o  2> $OUT.log.varann.e
if [[ $? == 0 ]]
then
        echo "Success GATK-ANNOTATE"
        #Annotate the variants further using ANNOVAR genomic annotation
        CMD="sh $BPATH/do_annovar_all.sh $OUT.annotatedvarnt.vcf "              ##gives output chr.vcf.annotatedvarnt.vcf.annovar.summary.genome_
        echo $CMD
        $CMD 1> $OUT.log.annovar.o  2> $OUT.log.annovar.e
        if [[ $? == 0 ]]
        then
                echo "Success ANNOVAR"
                        chkComplete=`grep "Success ANNOVAR" $Dir/LOG/calling.$sampleName.*.o | wc -l `
                        if [[ $chkComplete == $contigs ]]  ##chr 1-22,x y (& mt in case of WGS)
                        then
                                echo "Success ANNOVAR" > $Dir/LOG/calling.$sampleName.all.o
                                CMD="qsub -o $Dir/LOG/post.all.o -e $Dir/LOG/post.all.e -N post.calling -l mem=3G,time=4:: $BPATH/samtools_callin
                                echo $CMD
                                $CMD
                                echo "FINAL OUTPUT: $Dir/$sampleName/ALL"
                        fi
                        echo "End $0: `date`"
                        doexit 0
        else
                echo "Failed ANNOVAR"
                echo "End $0: `date`"
                doexit 1
        fi
else
        echo "Failed GATK-ANNOTATE"
        echo "End $0: `date`"
        doexit 1
fi