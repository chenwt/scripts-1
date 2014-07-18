#!/usr/bin/python
#J.HE
'''
Given a list of barcode, extract all the mutations in them
input: -b <barcode>
        -m <mutation maf file>
output: -o <mutations with sample names>
status: complete
'''

import sys,getopt
from collections import defaultdict
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -b <barcode list> \
        -m <mutation maf file from tcga>\
        -o <output>'
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hb:m:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-b"):
        barcodefile = arg
    elif opt in ("-m"):
        maffile = arg
    elif opt in ("-o"):
        output = arg

barcodeL = []
with(open(barcodefile)) as f:
    line = f.readline()
    while line:
        barcodeL.append(line.strip())
        line = f.readline()

# print len(barcodeL), barcodeL[:5]
infoIndex = [0, 4, 5, 6, 7, 15] 

outDict = defaultdict(list)
with(open(maffile)) as f:
    line = f.readline()
    while line:
        if line[0] == "#":
            line = f.readline()
            continue
        info = [a for id, a in enumerate(line.strip().split("\t")) if id in \
                infoIndex]
        if info[-1] in barcodeL:
            outDict["\t".join(info[1:5] + [info[0]] ) ]\
                    .append(info[-1]) 
        line = f.readline()

outf = open(output, 'w')

for k,v in sorted(outDict.items()):
    outf.write(k + "\t" + ";".join(v) + "\n")

outf.close()

print "[---END----]"
