#!/bin/urs/python
#input: <file with two column: TCGA barcode TCGA BAM FILE NAME>
#output: the file listed in the second column of the input file's name will be changes
#description: rename all downloaded TCGA bam file in current folder with tcga barcode<PID-SAMPLE-CODE> 


import os
#mapf = "summary_BRCA.txt";
mapf = "input_2_rename_bam.txt";
# mapf = "testfile.txt"
with open(mapf) as f:
	for line in f:
		[new,old] = line.split("\t");
		# print new + ":" + old;
		cmd = "mv " + old.rstrip() +" " + new.rstrip()
		# cmd = "mv " + new.rstrip() + " " + old.rstrip()
		print cmd
		os.system(cmd)
# print " ".join(new);
# print " ".join(old);

		# os.rename(filename, filename[7:]);


