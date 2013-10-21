#!/bin/bash

die () {
  echo -e >&2 "$@"
  exit 1
}

[ "$#" -ge 5 ] || die "BAFPLOT.sh requires five arguments:\n1) a pathname for a set of bed files, e.g., ../dbSNP135commonchrCHR.bed\nThe bed files at /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/DBSNPCommonSNPsbedfiles/ contain all SNPs from dbSNP with a frequency greater than 10%.\n2) a pathname for a set of bam files, e.g., ../CHR.bam.sorted.bam.noDup.bam.baq.bam\n3) a directory in which to save intermediate files and the BAF plot for the positions specified in the bed files in the specified bam files\n4) a title to place on the BAF plot (use underscores instead of spaces in the title, e.g., My_BAF_Plot)\n5) the file name for the plot.\nCHR in arguments 1 and 2 above will be replaced with 1..22.\nThe BAF plot will be saved in the specified directory.\nSample command:\n./BAFPLOT.sh /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/DBSNPCommonSNPsbedfiles/dbSNP135commonchrCHR.bed /ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/121016_ANDREA_GABRIELLE_12_HUMAN_GENOME_500M_PE100_HISEQ/PILOT_SAMPLES/AC3/15X/CHR.bam.sorted.bam.noDup.bam.baq.bam /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/AC3BAF AC3_BAF_PLOT AC3bafplot.pdf" 

bedfilemodel=$1
bamfilemodel=$2
workfolder=$3
plottitle=$4
plotname=$5

pileupfilemodel=$workfolder/chrCHR.pileup
freqfilemodel=$workfolder/chrCHR.freq
bafplot=$workfolder/$plotname
routputfile=${workfolder}/freqtobafplot.Rout


#use this version of R, or another version if it has libraries ggplot2, grid 
#and scaled installed
#alias R='/ifs/scratch/c2b2/ngs_lab/db2175/R-2.15.1/bin/R'

export BEDTOPILEUPSH='/ifs/scratch/c2b2/ys_lab/jh3283/net/BAFpipeline/bed_to_pileup.sh'
export MAKEPILEUP='/ifs/scratch/c2b2/ys_lab/jh3283/net/BAFpipeline/make_pileup.sh'

export PILEUPTOFREQSH='/ifs/scratch/c2b2/ys_lab/jh3283/net/BAFpipeline/pileup_to_freq.sh'
export PILEUPTOFREQPY='/ifs/scratch/c2b2/ys_lab/jh3283/net/BAFpipeline/pileup_to_freq.py'

export FREQTOBAFPLOTR='/ifs/scratch/c2b2/ys_lab/jh3283/net/BAFpipeline/freq_to_bafplot.R'

$BEDTOPILEUPSH $workfolder $bedfilemodel $bamfilemodel $pileupfilemodel
$PILEUPTOFREQSH $pileupfilemodel $freqfilemodel
/nfs/apps/R/2.15.1/bin/R CMD BATCH --slave "--args $freqfilemodel $bafplot $plottitle $routputfile" $FREQTOBAFPLOTR 
if [ $? -ne 0 ]; then
  echo "R did not exit properly.  Inspect $routputfile.  "
else 
  echo "Execution completed"
fi
