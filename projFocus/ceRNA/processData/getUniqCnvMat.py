#!/usr/bin/python
##J.HE
##Desp.: given CNV matrix, get uniq gene CNV by getting the larger one, input
#    should be sorted

import os,sys,getopt

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
print('Log file:'+ outlog)

def updateCNV(cnvlist1, cnvlist2):
    ''' update CNV value, take the one with bigger CNV"
    if len(cnvlist1) == len(cnvlist2):
        cntCnv = len(cnvlist1)
        cnvNew = [] * cntCnv 
        for i in range(cntCnv):
            if abs(cnvlist1[i]) > abs(cnvlist2[i]) :
                cnvNew.append(cnvlist1[i])
            else:
                cnvNew.append(cnvlist2[i])
    else:
        print "Warning: sample number not the same!"
    return cnvNew

genePre = ''
cnvsPre = 0
outputH = open(output, 'w')
headerFlag = 0
with open(input) as f:
    for line in f.readlines():
        geneNew, cnvsNew = line.strip().split("\t",1)
        if geneNew == genePre:
            cnvsNewlist = map(float,cnvsNew.split("\t"))
            cnvsPre = updateCNV(cnvsPre, cnvsNewlist)
        else:
            if geneNew == "barcode":
                genePre = geneNew
                cnvsPre = cnvsNew.split("\t")
                continue
            else:
                outputH.write(genePre + "\t" + \
                              "\t".join(map(str,cnvsPre)) + "\n") 
            genePre = geneNew
            cnvsPre = map(float, cnvsNew.split("\t"))

outputH.close()
print "SUCESS"
