#!/usr/bin/python
#J.HE
'''
Get reference genome coordinates length from NCBI36 hg
'''

t0 = time()
file = "/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/wu_build36.fasta"

hg_ncbi36_coord = []
with (open(file,'r')) as f:
        line = f.readline()
        while line:
            if line[0] == ">":
                hg_ncbi36_coord.append(line.strip().split(":")[-4:-1])
                line = f.readline()

t1 = time()
print t1-t0

fout = open( file + ".chrom_length", 'w')
for i in hg_ncbi36_coord:
    fout.write("\t".join(i) + "\n")
fout.close()

print "[-----END at " + str(time()) + "----]"
