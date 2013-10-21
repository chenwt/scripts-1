#!/bin/bash

die () {
  echo -e >&2 "$@"
  exit 1
}

[ "$#" -ge 2 ] || die "pileup_to_freq.sh requires two arguments:\n1) a pathname for a set of pileup files, e.g., ../AC3dbSNP135commonchrCHR.pileup\n2) an output pathname, e.g., ../AC3dbSNP135commonchrCHR.freq\nIn both cases CHR will be replaced with 1..22.\npileup_to_freq.sh calls pileup_to_freq.py. pileup_to_freq.py expects the pileup files to have no header and 6 columns: the chromosome, the position, the reference base, the number of reads, the pileup and the quality string.  pileup_to_freq.py generates a file with a 5 line header and then a table with 4 columns (starting with the column names: Chr, Pos, Totalreads, Altreads).\nSample command:\n./pileup_to_freq.sh /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/AC3DBSNPCommonSNPspileups/AC3dbSNP135commmonchrCHR.pileup /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/AC3DBSNPCommonSNPsfreqs/AC3dbSNP135commonchrCHR.freq"
[ $PILEUPTOFREQPY ] || die "Variable PILEUPTOFREQPY must be set to pathname of pileup_to_freq.py script.  Try /ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/pileup_to_freq.py"

pileupfilemodel=$1
outputfilemodel=$2

for i in {1..22}; do
  newpileupfile=${pileupfilemodel/CHR/$i}
  newoutputfile=${outputfilemodel/CHR/$i}
  if [ ! -f $newoutputfile ]; then
    echo "Writing frequency file for chromosome "$i
    python -B $PILEUPTOFREQPY $newpileupfile $newoutputfile
    if [ $? -ne 0 ]; then
      echo "python did not exit properly"
      mv $newoutputfile $newoutputfile.fail
    fi
  else
    echo "$newoutputfile already exists.  Not overwriting."
  fi
done
