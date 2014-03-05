#!/usr/bin/env python
# J.HE
# input: 1 MRs gene;
#		 a file has mutated gene names, one per line
#		 PrePPI database with genename as interactor.
# output: MRs, first order interacting proteins, mutated gene
#Desp:given a list of uniport gene, return it's linked proteins in preppi\
        # database, probability of 0.9

import os,getopt,sys
argv = sys.argv[1:]
input = ''
order = 1
usage = 'python " + sys.argv[0] + " -i <input>  '
example = 'python " + sys.argv[0] + " -i <input>'
preppi='/ifs/scratch/c2b2/ac_lab/jh3283/database/preppi/preppi_int_3col_genename.txt_600_90' 

try:
    opts,args = getopt.getopt(argv,"hi:l:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i"):
        input = arg
    elif opt in ("-l"):
        order = int(arg)


output = input + ".preppi"
print ('Script path'+ sys.argv[0])
print('Input file:' + input)
print('Output file:' + output)
print('order:',  order)

# get mutGene correlated protein
def intersect(a, b): #function to get intersection of two lists
	if ((len(a) > 0) & (len(b) > 0)):
		res = list(set(a) & set(b))
	else:
		res=[]
	return res
gene = []
with open(input) as f:
    for line in f.readlines():
        gene.append(line.strip().split()[0])
print "input gene" + "\t" + str(len(gene))

def getInteractor(gene, network):
    ppi = []
    cnt = 0 
    with open(network) as f:
    	for line in f.readlines():
            [p1, p2, val] = line.split("\t",2)[:3]
            if p1 in gene or p2 in gene and not [p1,p2] in ppi: 
                ppi.append([p1,p2])
                cnt = cnt + 1
    print  "interaction\t" + str(cnt)
    return ppi

ppi = getInteractor(gene,preppi)

if order == 1 :
    outputH = open(output,'w')
    for [p1,p2] in ppi:
        outputH.write( p1 + "\t" + p2 + "\n")
    outputH.close()

elif order == 2:
    gene = [p for pp in ppi for p in pp ] 
    outputH = open(output + "." + str(order) + ".gene", 'w')
    outputH.write("\n".join(gene))
    outputH.close()
    ppi = getInteractor(gene,preppi)
    outputH = open(output,'w')
    for [p1,p2] in ppi :
        outputH.write( p1 + "\t" + p2 + "\n")
    outputH.close()
else:
    print "underdevelopment! "

print "SUCESS"
