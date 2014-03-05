#!/usr/bin/python
#J.HE
##get matrix using  output form getMAF.py
## use the same genelist<or smaller subset> one input
## output granularity base on the gene


##-----------------funcS
def matchSample(smp,smplist):
    """given a string,check if it matches the sample, and return the sample it
    matches,blur version"""
    if smp in smplist:
        return smp
    else:
        return False

class MutGroup(chr,ps,pe,gname,maf,barcode):
    """docstring for Mutation
    store somatic mutation
    """
    def __init__(self, chr, ps, pe, gname, maf, barcode):
        self.chr = chr
        self.ps  = long(ps) # starting point 
        self.pe  = long(pe) # ending point, mid position of  mutations sites in one mut cloud
        self.gname  = gname # gene name  
        self.maf = float(maf) # largest maf in this mutation cloud
        self.barcode = barcode
        tempMutCode = ":".join(map(str, [mut2.chr, mut2.ps, mut2.pe]))
        self.mutCode.append([tempMutCode])

    def updateMut(self, mut2):
        """ check mut2 in tied to self by distance, REG = 100 bp
        return updated mut1, """
        if self.gname == mut2.gname and self.barcode = mut2.barcode:
            if  self.ps > mut2.pe + REG or self.pe < mut2.ps + REG:
                return 0 
            else: 
                self.ps = min( self.ps, mut2.ps)
                self.pe = max( self.pe, mut2.pe)
                self.maf = max( self.maf, mut2.maf)
                tempMutCode = ":".join(map(str, [mut2.chr, mut2.ps, mut2.pe]))
                self.mutCode.append([tempMutCode])
                return 1 
        else: 
            return 0
    def inGene(self, gene):
        if self.gname == gene:
            result 1
        else:
            result 0

##-----------------funcE

import numpy as np
import os,getopt,sys

REG = 100.0
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " \
        -i <input>  \
        -g <genelist,first colum genename> \
        -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hi:g:o:",["ifile=","genelist=","ofile="])
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
        glistFile = arg
    elif opt in ("-s","--samplelist"):
        slistFile = arg
    elif opt in ("-o","--ofile"):
        output = arg
        outlog = output + ".log"
print ('Script path'+ sys.argv[0])
print('Input file:' + input)
print('Output file:'+ output)
print('Log file:'+ outlog)

# geneName        mutCode MAF     barcode
# INPP5A  10:133730291:133730292  0.292682926829  H_LS-E2-A15H-01A-11D-A19H-09
#load genename
glist = []
with open(glistFile) as f:
    for line in f.readlines():
        glist.append(line.strip().split()[0])
cntG = len(glist)
print "number of input genes\t" + str(cntG)

#load sampleid
slist = []
with open(slistFile) as f:
    for line in f.readlines():
        slist.append(line.strip().split()[0])
cntS = len(slistFile)
print "number of input samples\t" + str(cntS)

#load maf file
cnt=0
mafArray = np.zeros( (cntG,cntS),dtype=float)

geneMutGrpDict = {} # 
resultDict = {} #dict of Mutation class, value is barcode list 
with open(input, buffering = 100000) as f:
    f.next()
    for line in f.readlines():
        tempGene, tempChr, tempPs, tempPe, tempMaf, tempBarode = line.strip().split("\t|:")
        if tempGene in glist and tempBarcode in slist:
            print "\t".join([ tempGene, tempChr, tempPs, tempPe, tempMaf, tempBarcode]) 
            tempMut = MutGroup(  tempChr, tempPs, tempPe, tempGene, tempMaf,\
                               tempBarcode)
            tempGeneMutlist = geneMutGrpDict.get(tempGene, 0)
            if tempGeneMutlist != 0 :
                for a in tempGeneMutlist: 
                    if not a.updateMut(tempMut) :
                        tempGeneMutlist.append(tempMut)
                    else:
                        break ## break once find one group                      
            else :
                mutGrpDict[tempGene] = [tempMut] 
        else:
            continue

outputH = open(output,'w+')
outputH.write("gene\tmutloc\t" + "\t".join(glist) + "\n")
for g in geneMutGrpDict.keys() :
    for m in geneMutGrpDict[g]:
        outputH.write(g + "\t" +  
