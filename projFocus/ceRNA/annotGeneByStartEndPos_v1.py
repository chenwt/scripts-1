#!/usr/bin/env python
#J.HE
#Desp:created for projFocus to annot each somatic mutated gene by gene starting and ending position.Dec 17,2013
#     this code need to be run before connect with expression data, find the promoter region somatic mutation gene
#input:<file:  genenlsit > 
#output:<file: genelist.startend> 
#TODO:

dbfile = "/ifs/scratch/c2b2/ac_lab/jh3283/database/projFocusRef/refset_gene_start_end.tsv.sorted.uniq_resortedCol_geneSingleStartEnd"
import os,getopt,sys
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i","--ifile"):
        input = arg
    elif opt in ("-o","--ofile"):
        output = arg
        outlog = output + ".log"
print ('Script path'+ sys.argv[0])
print('Input file:' + input)
print('Output file:'+ output)

gene = []
with open(input) as f:
    for line in f.readlines():
        gene.append(line.strip())
print "input gene" + str(len(gene))

outputH = open(output,'w')
with open(dbfile) as f:
    for line in f.readlines():
        tempGene, vars = line.strip().split("\t",1)
        if len(gene) > 0  :
            if tempGene in gene:
                outputH.write(tempGene + "\t" + vars + "\n") 
                gene.remove(tempGene)
        else:
            break
outputH.close()
print "SUCESS"
