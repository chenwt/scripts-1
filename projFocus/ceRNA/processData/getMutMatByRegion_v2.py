#!/usr/bin/python
#J.HE
## get matrix using  output form splitByKey.py
## use the same genelist<or smaller subset> one input
## output 1. mat of group by barcode maf, 2 mapping file group, mutation, 
#      barcode base on the gene
import os,getopt,sys
import re 
import operator
from subprocess import Popen, PIPE

usage = 'python  + sys.argv[0] + \n \
        -h usage and  example  \n\
        -i --idir inputfile dirname, the output of splitByKey.py \n\
        -g --genelist" genelistFile, first column must be the genename\n\
        -s --barcodelist barcodeFile, exact match as in the input file\n\
        -d--distance <optional> DISTANCE for grouping, defalut 100 bp \n\
        -o --ofile  output file name\n'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'

##-----------------funcS
class Mut():
    def __init__(self, gene,chr, ps, pe,maf,barcode):
        self.chr = chr
        self.ps = long(ps)
        self.pe = long(pe)
        self.maf = maf
        self.gene = gene
        self.barcode = barcode
    def show(self):
        print "\t".join(map(str, \
        [self.chr,self.ps,self.pe,self.maf, self.barcode]))
  
def initOutValue(mut,Order):
    outValue = [0] *  len(Order) 
    outValue[:4] = [mut.gene, mut.chr, mut.ps, mut.pe]
    tempIndx = Order.index(mut.barcode)
    outValue[tempIndx] = max(outValue[tempIndx], mut.maf)
    return outValue
    
def initOutList(mut):
    outList = [0] * 5
    outList[:4] = [mut.gene, mut.chr, mut.ps, mut.pe]
    outList[4] = ['0+'+ str(mut.pe-mut.ps)]
    return outList
##-----------------funcE

ERR = "#ERROR at "
SUC = "#Sucess of getMutMatByRegion "   
DISTANCE = 100 
argv = sys.argv[1:]

try:
    opts,args = getopt.getopt(argv,"hi:d:s:o:",\
    ["ifile=","distance=","barcodelist=","ofile="])
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i","--ifile"):
        tempDir   = arg
        tempDir   = os.path.abspath(tempDir)
    elif opt in ("-s","--barcodelist"): 
        barcodeFile = arg
        barcodeFile = os.path.abspath(barcodeFile)
    elif opt in ("-d", "--distance"):
        DISTANCE = int(arg)
    elif opt in ("-o","--ofile"):
        output = arg ###abs path
        output2 = output+ ".grp_muts_Mapping"
print('Script path'+ sys.argv[0])
print('Input file:' + tempDir )
print('Output file:'+ output)

##---getting barcode information
with open(barcodeFile) as f:
    barcodeList = list([l.strip() for l in f.readlines()])

##---getting genefile information
try:
    subfiles = os.listdir(tempDir)
except:
    print "Error getting subfiles.."
    sys.exit()

##---output file
outOrder = ["geneName","mutChr","mutRegStart","mutRegEnd"]+ barcodeList
outputH = open(output, 'w+')
outputH.write("\t".join(outOrder) + "\n")
outputH2 = open(output2, 'w+')
outputH2.write("\t".join(["geneName","mutRegChr","mutRegStart","mutRegEnd",\
                "mut1:mut2:..."]) + "\n")
    
for gfile in subfiles:
    with open(tempDir + "/" + gfile) as f:
        mutList = []
        for line in f.readlines():
            v1,v2,v3,v4,v5,v6 =  re.split(r"\t|:",line.strip())
            tempMut = Mut(v1,v2,v3,v4,v5,v6)
            mutList.extend([tempMut])
    # sort and output   
    mutLen = len(mutList)
    if mutLen > 1:
        mutPrev = mutList[0]
        mutList.sort(key = lambda x: x.ps, reverse=False)
        outValue = initOutValue(mutPrev,outOrder)
        outList  = initOutList(mutPrev)
        for mutCrt in mutList[1:]: 
            if mutCrt.ps - mutPrev.pe > DISTANCE:
                #print "\t".join(map(str, outValue))
                outputH.write("\t".join(map(str, outValue)) + "\n")
                outList[4] = ":".join(outList[4])
                outputH2.write("\t".join(map(str, outList)) + "\n")

                mutPrev = mutCrt
                outValue = initOutValue(mutPrev,outOrder)     
                outList  = initOutList(mutPrev)
            else:
                outValue[3] = max(mutCrt.pe,mutPrev.pe) 
                outValue[2] = min(mutCrt.ps,mutPrev.ps)                               
                tempIndx = outOrder.index(mutCrt.barcode)
                outValue[tempIndx] = max(outValue[tempIndx], mutCrt.maf)
                temp = str(mutCrt.ps - outValue[2]) + "+" +str(mutCrt.pe - mutCrt.ps)
                outList[4] = outList[4] + [temp] 
                mutPrev = mutCrt
                
    elif mutLen == 1:
        mutPrev = mutList[0]
        outValue = initOutValue(mutPrev,outOrder)
        outList  = initOutList(mutPrev)
        outputH.write("\t".join(map(str, outValue)) + "\n")     
        outList[4] = ":".join(outList[4])
        outputH2.write("\t".join(map(str, outList)) + "\n")
    else:
        print "no mutation for gene" + gfile

##---cleaning up
# cmd = "rm -rf "+ tempDir
# p = Popen(cmd, stdout=PIPE,stderr=PIPE, shell=True)
# p.communicate()

print SUC
