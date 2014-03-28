#!/usr/bin/env python
#J.HE
#Desp: used to get single start and end for each gene instead of transcripts for entrez_annotation_hg19.txt file
#input: entrez_annotation_hg19.txt
#output: output also tss_tse 
##COMMENT: inputfile well sorted by chromosome

import os,re
import sys,getopt
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
        inpFile = arg
    elif opt in ("-o","--ofile"):
        outFile = arg
print('Script path'+ sys.argv[0])


def getBigCnv(vls):
    if abs(vls[0]) >= abs(vls[1]):
        return vls[0]
    else: 
        return vls[1]


outFileHandler = open(outFile,'w')
cntLine = 0 
geneInfoOld = []
with open(inpFile) as f:
    line = f.readline()
    while line:
        if re.match("barcode",line) or re.match('^gen',line) :
            geneInfoOld = line.strip().split("\t", 1)
        else:
            [geneId, val]= line.strip().split("\t", 1)
            geneInfoCrt = [geneId,val]
            if not geneInfoOld: 
                geneInfoOld = geneInfoCrt 
            elif geneInfoOld and geneInfoOld[0] == geneId : 
                oldvalTemp = map(float, geneInfoOld[1].split("\t"))
                crtvalTemp = map(float, val.split("\t"))
                # print map(getBigCnv, zip(oldvalTemp, crtvalTemp) )
                try:
                    geneInfoOld[1] = "\t".join(map(str,\
                                map(getBigCnv, zip(oldvalTemp, crtvalTemp) ))) 
                except:
                    print "ERROR in length of value"
                    sys.exit()
            elif geneInfoOld and not geneInfoOld[0] == geneId : 
                outFileHandler.write( "\t".join(geneInfoOld) + "\n" )
                geneInfoOld = geneInfoCrt 
            else:
                print "ERROR"
                sys.exit()
        cntLine = cntLine + 1
        line = f.readline()
    print "line processed:\t" + str(cntLine)

outFileHandler.close()
print "#----END---"

