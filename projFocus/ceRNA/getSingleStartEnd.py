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

outFileHandler = open(outFile,'w')
cntLine = 0 
geneInfoOld = ["gene","chrom","start","end","strand"]
with open(inpFile) as f:
    line = f.readline()
    while line:
        if re.match("^#",line) or re.match('^gen',line) :
            pass
        else:
            [geneId, chrom, start, end, strand] = line.strip().split("\t",4)
            geneInfoCrt = [geneId, chrom, start, end, strand]
            if not geneInfoOld: 
                geneInfoOld = geneInfoCrt 
            elif geneInfoOld and geneInfoOld[0] == geneId : 
                if geneInfoCrt[4] == '-':
                    if long(geneInfoCrt[2]) > long(geneInfoOld[2]):
                        geneInfoOld[2]  = start
                        geneInfoOld[3]  = end
                    else:
                        pass
                elif geneInfoCrt[4] == '+' :
                    if long(geneInfoCrt[2]) < long(geneInfoOld[2]) :
                        geneInfoOld[2]  = start
                        geneInfoOld[3]  = end
                    else:
                        pass
                else:
                    pass 
            elif geneInfoOld and not geneInfoOld[0] == geneId : 
                outFileHandler.write( "\t".join(geneInfoOld) + "\n" )
                geneInfoOld = geneInfoCrt 
        cntLine = cntLine + 1
        line = f.readline()
    print "line processed:\t" + str(cntLine)

outFileHandler.close()
print "#----END---"

