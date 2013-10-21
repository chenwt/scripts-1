#!/bin/bash

die () {
  echo -e >&2 "$@"
  exit 1
}

[ "$#" -ge 4 ] || die "bed_to_pileup.sh requires 4 arguments:\n1) the working folder (for temp files etc...)\n2) a pathname for a set of bed files, e.g., ../dbSNP135commonchrCHR.bed\n3) a pathname for a set of bam files, e.g., ../CHR.bam.sorted.bam.noDup.bam.baq.bam\n4) an output pathname, e.g., ../CHR.pileup\nIn all cases CHR will be replaced with 1..22.\nbed_to_pileup.sh calls make_pileup.sh, which calls samtools to create pileup files for the locations specified in the bed files from the bam files.\nSample command:\n./bed_to_pileup.sh /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/DBSNPCommonSNPsbedfiles/dbSNP135commonchrCHR.bed /ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/121016_ANDREA_GABRIELLE_12_HUMAN_GENOME_500M_PE100_HISEQ/PILOT_SAMPLES/AC3/15X/CHR.bam.sorted.bam.noDup.bam.baq.bam /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/AC3DBSNPCommonSNPspileups/AC3dbSNP135commmonchrCHR.pileup"
[ $MAKEPILEUP ] || die "Variable MAKEPILEUP must be set to pathname of make_pileup.sh script.  Try /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/make_pileup.sh."  
[ $REF ] || die "Variable REF must be set to pathname of reference."  

workfolder=$1
bedfilemodel=$2
bamfilemodel=$3
outputfilemodel=$4

for i in {1..22}; do
  newbedfile=${bedfilemodel/CHR/$i}
  newbamfile=${bamfilemodel/CHR/$i}
  newoutputfile=${outputfilemodel/CHR/$i}
  if [ ! -f $newoutputfile ]; then
    echo "Writing pileup file for chromosome "$i "using reference file "$REF
    $MAKEPILEUP $workfolder $newbamfile $REF $newbedfile $newoutputfile 
    if [ $? -ne 0 ]; then
      echo "make_pileup.sh did not exit properly"
      mv $newoutputfile $newoutputfile.fail
    fi
  else
    echo "$newoutputfile already exists.  Not overwriting."
  fi
done
