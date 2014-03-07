#!/usr/bin/python
#J.HE
##get matrix using  output form getMAF.py
## use the same genelist<or smaller subset> one input
## output 1. mat of group by barcode maf, 2 mapping file group, mutation, 
#      barcode base on the gene
import os,getopt,sys
import re 
import operator
from subprocess import Popen, PIPE

usage = 'python  + sys.argv[0] + \n \
        -h usage and  example  \n\
        -i --ifile inputfile name, the output of step-6_getMutMaf.py \n\
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
SUC = "#Sucess "   
DISTANCE = 100 
argv = sys.argv[1:]
# ##----testS
# input = 'test.mat'
# genelistFile = "gene.list"
# output = 'test.matrix'
# output2 = output+ ".grp_muts_Mapping "
# barcodeFile = "barcode.list"
# ####----testE

try:
    opts,args = getopt.getopt(argv,"hi:g:d:s:o:",\
    ["ifile=","genelist=","distance=","barcodelist=","ofile="])
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i","--ifile"):
        input = arg
    elif opt in ("-g","--genelist"):
        genelistFile = arg
    elif opt in ("-s","--barcodelist"): 
        barcodeFile = arg
    elif opt in ("-d", "--distance"):
        DISTANCE = int(arg)
    elif opt in ("-o","--ofile"):
        output = arg
        output2 = output+ ".grp_muts_Mapping"
print('Script path'+ sys.argv[0])
print('Input file:' + input)
print('Output file:'+ output)

##split file by input genelistfile
tempDir =  output+ "-tempDir"
# cmd = "~/tools/python/Python_current/python " + "~/bin/splitByKey.py" +\
cmd = "python " + "~/bin/splitByKey.py" +\
        " -i " + input + " -s 1000000 " + " -k " + genelistFile \
        + " -o " + tempDir
print cmd
p1 = Popen(cmd, stdout=PIPE,stderr=PIPE, shell=True)
err = p1.communicate()[1]
if err:
    print ERR + "split file by genes" 
    sys.exit() 

##---getting barcode information
with open(barcodeFile) as f:
    barcodeList = list([l.strip() for l in f.readlines()])

CWD = os.path.abspath(".")
try:
    subfiles = os.listdir(CWD + "/" + tempDir)
except:
    print "Error getting subfile.."
    sys.exit()

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
cmd = "rm -rf "+ tempDir
p = Popen(cmd, stdout=PIPE,stderr=PIPE, shell=True)
p.communicate()

print SUC
