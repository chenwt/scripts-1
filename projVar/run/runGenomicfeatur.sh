#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

function installHgwiggle {
  curl -O http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
  mkdir -p ~/bin/${MACHTYPE}
}

function download {
  curl -O http://hgdownload.cse.ucsc.edu/gbdb/hg19/multiz46way/phastCons46wayPrimates.wib
  
  curil -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/phastCons46wayPrimates.txt.gz
}


###---replication time


