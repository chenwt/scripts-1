#!/usr/bin/python
#J.HE
#Desp.: given cnv matrix, meth matrix, som matrix,sample names,gene list return the gene sample
# namelist ( gene\tsample1;sample2;sample3)


import re
import sys,getopt
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " -b <barcodelist> -g <genelist> -c <cnvmat> \
-m <methmat> -s <sommat>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'

try:
    opts,args = getopt.getopt(argv,"hb:g:c:m:s:o:",["ifile=","ofile="])
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-b","--ifile"):
        fbcode = arg
    elif opt in ("-c","--ifile"):
        fcnv = arg
    elif opt in ("-m","--ifile"):
        fmeth  = arg
    elif opt in ("-s","--ifile"):
        fsom = arg
    elif opt in ("-o","--ofile"):
        output = arg
        outlog = output + ".log"
print('Script path:\t"'+ sys.argv[0])
print('Input file:\t' + input)

bcodelist = [] 
with open(fbcode) as f:
    line = f.readline()
    while line:
        bcodelist.append(re.sub("-", ".", line.strip()))
        line = f.readline()

gintlist = {}
cnt = 0
with open(fcnv) as f:
    line = f.readline()
    while line:
        if re.match("^gene|barcode",line):
            bcCnvlist = re.sub("-",".", line.strip()).split("\t")
        else: 
            # tempIdxCnvfree = filter(lambda x: float(x) != 0.0,\
                                      # line.strip().split("\t",1)
            templine = line.strip().split("\t")
            tempIdxCnvfree = [ idx + 1 for (idx, val) in \
                              enumerate(templine[1:]) if \
                               float(val) == 0.0 ]  
            cnt = cnt + 1
            tempGintBcode = map(bcCnvlist.__getitem__, tempIdxCnvfree)
            # if gintlist.get(templine[0],0):
            gintlist[templine[0]] = tempGintBcode 
        line = f.readline()
print "cnt:\t" + str(cnt)
cnt = 0 
with open(fmeth) as f:
    line = f.readline()
    while line:
        if re.match("^gene|barcode",line) or re.findall("\.01A", line):
            bcMethlist = re.sub("-",".", line.strip()).split("\t")
        else: 
            templine = line.strip().split("\t")
            tempIdxMethfree = [ idx for (idx, val) in \
                               enumerate(templine[1:]) if \
                               int(val) == 0 ]  
            tempGintBcode = map(bcMethlist.__getitem__, tempIdxMethfree)
            if gintlist.get(templine[0],0):
                gintlist[templine[0]] = [ ss for ss in tempGintBcode if ss in \
                                      gintlist.get(templine[0],0) ] 
        cnt = cnt + 1
        line = f.readline()
print "meth:\t" + str(cnt)

outputH = open(output + "_CnvMethFree" , 'w')
outHeader = "geneName" + "\t" + "gintSamples" + "\n" 
outputH.write(outHeader) 
for key, val in gintlist.items():
    outputH.write(key + "\t" + ";".join(val) + "\n") 
outputH.close()

cnt = 0 
with open(fsom) as f:
    line = f.readline()
    while line:
        if re.match("^gene|barcode",line) or re.findall(".01A", line):
            bcSomlist = re.sub("-",".", line.strip()).split("\t")
            bcSomlist = map(lambda x:x[6:17], bcSomlist[1:]) 
        else: 
            templine = line.strip().split("\t")
            # print templine 
            tempIdxSomfree = [ idx  for (idx, val) in \
                               enumerate(templine[1:]) if \
                               int(val) == 0 ]  
            tempGintBcode = map(bcSomlist.__getitem__, tempIdxSomfree)
            if  gintlist.get(templine[0],0):
                gintlist[templine[0]] = [ ss for ss in tempGintBcode if ss in \
                                      gintlist.get(templine[0],0) ] 
        cnt = cnt + 1
        line = f.readline()
print "som:\t" + str(cnt)
outputH = open(output + "_CnvMethSomFree" , 'w')
outHeader = "geneName" + "\t" + "gintSamples" + "\n" 
outputH.write(outHeader) 

for key, val in gintlist.items():
    outputH.write(key + "\t" + ";".join(val) + "\n") 
outputH.close()


