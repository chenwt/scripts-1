#!/usr/bin/python
#J.HE
'''
Desp.: generate genome wide random mutation,
input: -n  number of mutation
        -s  positive mutation file,(will be excluded from randome)
output: randome file pickle file
status : complete
'''

from collections import defaultdict
import sys,getopt
import random

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -i <positive mutation file>  \
        -n <number of random mutation> -o <output>'

try:
    opts,args = getopt.getopt(argv,"hs:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n"  
        sys.exit()
    elif opt in ("-s"):
        somfile = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + somfile)
print('Output file:\t'+ output)

# functions 

def loadPosMut(file):
    '''
    load positive mutations with first colum as chrom, second as position
    '''
    posMutD = defaultdict(list)
    chromL = map(str, range(1,23) ) + ['X', 'Y', 'x', 'y']
    cnt = 0 
    with(open(file)) as f:
        line = f.readline()
        while line:
            info1, info2 = line.strip().split()[:2]
            if info1.replace("chr","") in chromL:
                posMutD[info1].append(info2)
                cnt = cnt + 1
            line = f.readline()
    return posMutD, cnt

def genRandPtMut(n, chromLenD):
    '''
    generate randome mutations giveing a number
    '''
    rdMut = defaultdict(list)
    for i in range(n):
        rdChr = random.randint(1,24) 
        rdPos = random.randint(1,chromLenD[rdChr])
        rdMut[rdChr].append(rdPos)
    return rdMut

def rmSomFromRandMut(rdMut, smMut):
    '''
    work on list of each algorithm
    '''
    for k,vL in rdMut.items():
        if smMut.get(k, 0) == 0: 
           continue 
        vNewL = [ v for v in vL if v not in smMut[k]]
        if vNewL:
            rdMut[k] = vNewL
        else:
            del rdMut[k]     
    return rdMut
#     return dict((k,v) for k,v in rdMut.iteritems() if v)


###------start 
# store chrom length
file = "/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/wu_build36.fasta.chrom_length"

chromL = map(str, range(1,23) ) + ['X', 'Y']
chromLenD = defaultdict(int)
with open(file) as f:
    line = f.readline()
    while line:
        if line[0] in ['M', 'N']:
            line = f.readline()
            continue
        chrom, ps, pe = line.strip().split("\t")
        chromLenD[chromL.index(chrom) + 1] = int(pe) + 1
        line = f.readline()


somMut1D , n = loadPosMut(somfile)
n = int( n * 1.2)

rdMut1D = genRandPtMut(n, chromLenD)

## remove somatic from random 
randMutD = rmSomFromRandMut(rdMut1D, somMut1D)

## output 

outf = open(output, 'w')

for k, vL in randMutD.items():
    # print 
    # if k == 23:
    #     k = 'X'
    # elif k == 24:
    #     k = 'Y'
    # else:
    #     k = str(k)
    for v in vL:
        v = str(v)
        outf.write( "chr" + chromL[k-1]  + "\t" + v + " " + v + "\n" )

outf.close()

# import pickle 
# pickle.dump(randMutD, open(output  +".pkl","wb"))

print "[----END----]"

