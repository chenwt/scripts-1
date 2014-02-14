#!/bin/bash
#By: J.He
#TODO: 
#$-cwd
  
 REGIONFILE=$1
 inputBAM=$2
 outputBAMDIR=$3
 sampleName=`echo $inputBAM|awk 'BEGIN{FS="/"}{print $NF}'`
 cd $outputBAMDIR
 if [ ! -d $outputBAMDIR/splitLog ]; then mkdir $outputBAMDIR/splitLog ; fi
 logDir=$outputBAMDIR/splitLog

 cnt=1
 while read line
 do
   echo -e "chrom\tstart\tend\tcomment" > $outputBAMDIR/temp-region_$cnt.bed
   echo -e $line >> $outputBAMDIR/temp-region_$cnt.bed
   echo -e '#!/bin/bash' > $outputBAMDIR/temp-region_$cnt.sh
   echo "samtools view -hb -L $outputBAMDIR/temp-region_$cnt.bed $inputBAM > $outputBAMDIR/split_$cnt.$sampleName " >> $outputBAMDIR/temp-region_$cnt.sh
   echo "samtools index $outputBAMDIR/split_$cnt.$sampleName.sorted" >>$outputBAMDIR/temp-region_$cnt.sh 
   subBAM=`readlink -f $outputBAMDIR/split_$cnt.$sampleName.sorted`
   echo "echo $subBAM > $outputBAMDIR/split_$cnt.$sampleName.list " >>$outputBAMDIR/temp-region_$cnt.sh
   chmod +x $outputBAMDIR/temp-region_$cnt.sh 
   qsub -l mem=4g,time=8:: -e $logDir -o $logDir -cwd -N split_$cnt $outputBAMDIR/temp-region_$cnt.sh >> $outputBAMDIR/qsubSplit.log
   tail -1 $outputBAMDIR/qsubSplit.log 
   let cnt=$cnt+1
 done < $REGIONFILE


