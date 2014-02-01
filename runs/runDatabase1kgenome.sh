#!/bin/bash
#By: J.He
#TODO: 


#curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
# curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.sites.vcf.gz 
# gzip -d ALL.2of4intersection.20100804.sites.vcf.gz

###get slice vcf 
#tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 2:39967768-39967768

exCol(){
  n=$1
  file=$2
  out=$3
  awk -v n=$n '!/^#/&&$n!="."{print $n}' $file >> $out

}

exCol 3 ALL.2of4intersection.20100804.sites.vcf  ALL.2of4intersection.20100804.sites.rsID.txt
