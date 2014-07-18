#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

##
## This will requiring creating a directory for the file path in your directory. For example, if you were working on the phastCons46way for the hg19 database you would take the following steps:
##  1. Add a line to .hg.conf to point to where you are working
##gbdbLoc1=/home/usrName/work/
##  2. Obtain the wib file
##rsync -aP rsync://hgdownload.cse.ucsc.edu/gbdb/hg19/multiz46way/phastCons46way.wib  .
## 3. From that directory you pointed to, create a directory for the wib file
## mkdir -p hg19/multiz46way/
## 4. Move the file to that location.
## mv phastCons46way.wib hg19/multiz46way/
## 5. You can now perform operations like
## hgWiggle -db=hg19 -chr=chrM phastCons46way
## hgWiggle -db=hg19 -bedFile=bedFile phastCons46way
##

bedFile=$1
cd /ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/genomicFeature/conservation 
# ~/tools/hgWiggle -db=hg19 -bedFile=$bedFile -doBed phastCons46wayPrimates > $2
~/tools/hgtools/hgWiggle -bedFile=$bedFile  phastCons46wayPrimates > $1.cons
echo "[END]"
