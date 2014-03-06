#!/usr/bin/python
#J.HE
##get matrix using  output form getMAF.py
## use the same genelist<or smaller subset> one input
## output 1. mat of group by barcode maf, 2 mapping file group, mutation, 
#      barcode base on the gene


##-----------------funcS
def matchSample(smp,smplist):
    """given a string,check if it matches the sample, and return the sample it
    matches,blur version"""
    if smp in smplist:
        return smp
    else:
        return False

class MutGroup(gname,chr,ps,pe):
    """docstring for Mutation
    store somatic mutation
    """
    def __init__(self,gname, chr, ps, pe):
        self.gname  = gname # gene name  
        self.chr = chr
        self.ps  = long(ps) # starting point 
        self.pe  = long(pe) # ending point, mid position of  mutations sites in one mut cloud
        # self.maf = float(maf) # largest maf in this mutation cloud
        # self.barcode = barcode
        # tempMutCode = ":".join(map(str, [mut2.chr, mut2.ps, mut2.pe]))
        # self.mutCode.append([tempMutCode])

    def updateMut(self, mut2):
        """ check mut2 in tied to self by distance, REG = 100 bp
        return updated mut1, """
        if self.gname == mut2.gname : 
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
    def hasMut(self, mut2):
        if self.chr == mut2.chr:
            if self.ps > mut2.pe or self.pe < mut2.ps:
                return 0 
            else:
                return 1
        else:
            return 0
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
barcodelist = []
# cntS = len(slistFile)
# print "number of input samples\t" + str(cntS)

#load maf file from get_MAF.py output file
#get new group
cnt=0
# mafArray = np.zeros( (cntG,cntS),dtype=float)

def renewLocus(l1,l2):
    l1c, l1s, l1e = l1.split(":")
    l2c, l2s, l2e = l2.split(":")
    if l1c == l2c:
        ls = min(long(l1s), long(l2s))
        le = max(long(l1e), long(l2e))
        return ":".join(map(str,[l1c,ls,le]))
    else:
        print "Error in compare locus, different chroms"
        return 0

def compareLocus(l1,l2):
    l1c, l1s, l1e = l1.split(":")
    l2c, l2s, l2e = l2.split(":")
    l1s,l1e,l2s,l2e  = map(long, [l1s,l1e,l2s,l2e])
    if l1c == l2c:
        if l1s > l2e or l1e < l2s:
            return 0
        else:
            return 1
    else:
        return 0

# geneMutGrpDict = {} # 
grpBarcodeDict = {} ##dict of 
grpMafDict = {} #dict of mutGroup, value is barcode list 
#require 20M memory at least
### first scan to create group
with open(input, buffering = 10000000) as f:
    f.next()
    for line in f.readlines():
        tempGene, tempChr, tempPs, tempPe, tempMaf, tempBarode = line.strip().split("\t|:")
        if tempGene in glist and tempBarcode in slist:
            print "\t".join([tempGene, tempChr, tempPs, tempPe, tempMaf, tempBarcode]) 
            tempGrpMut = MutGroup(tempGene, tempBarcode)
            crtLocus = ":".join([tempChr,tempPs,tempPe])

            print "debug", crtLocus
            if not tempBarcode in barcodelist:
                ##update barcode list
                barcodelist.append(tempBarcode)
            tempGrpMaflist = grpMafDict.get(tempGrpMut, 0)
            preLocus = tempGrpMaflist[0]

            print "debug", preLocus
            if tempGrpMaflist != 0 : ##  find a existing [gene, barcode] bundle
                flagGrouped = 0 
                for a in tempGrpMaflist[1:]: 
                    if compareLocus(preLocus,crtLocus):
                        grpMafDict[tempGrpMut][0] = renew(preLocus,crtLocus)
                        grpMafDict[tempGrpMut].append(tempMaf)
                        flagGrouped = 1
                        ##need to update grpBarcodeDict
                    else:
                        continue ## skip to next once find one group                      
                if flagGrouped == 0:
                    grpMafDict[tempGrpMut][0] =
            else :
                grpMafDict[[tempGene,tempBarcode]] = tempMut 
        else:
            continue

outputH = open(output,'w+')
outputH.write("gene\tmutloc\t" + "\t".join(glist) + "\n")
for g in geneMutGrpDict.keys() :
    for m in geneMutGrpDict[g]:
        outputH.write(g + "\t" +  
