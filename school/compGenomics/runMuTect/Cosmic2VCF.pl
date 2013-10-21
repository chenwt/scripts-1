#!\usr\bin\perl -w
#Author: Jing He
#Date: Mar.28, 2013
#Last Update: Mar 29, 2013
#Usage:
#Example:

## use this to get corresponding columns in cosmic tsv file
# cat CosmicMutantExport_v54_080711.tsv | cut -f12,18 > pre_CosmicMutantExport_v54_080711.tsv 
# tail -n +2 /ifs/scratch/c2b2/ys_lab/jh3283/ref/test/pre_CosmicMutantExport_v54_080711.tsv|head -100 > /ifs/scratch/c2b2/ys_lab/jh3283/ref/test/test.tsv
use warnings;
# use strict;

#-----------Declaration---------
my $inFile ="/ifs/scratch/c2b2/ys_lab/jh3283/ref/test/test.tsv";
my $ouFile ="/ifs/scratch/c2b2/ys_lab/jh3283/ref/test/test.vcf";
my $line ="";
my @sline = ();
my @subString=();

open INF,$inFile;

while(<INF>) {
	$line =$_;
	@sline = split(/\W+/,$line);
	 # = ~/\d+([a-zA-Z]+)$/; 
	print "\n";
	for (my $i =0; $i < $#{$sline[0]} ;$i++){
		# print @{$sline[$i];
		if($sline[$i][1] =~ /\d+([a-zA-Z]+)$/){
			print "$1\t";
		}
	}
}
print "\n";
close(INF);

