#!/usr/bin/python
#J.HE
#Desp.: generate counts for venn diagram for each patients, all indels will be
# categorized into 7 groups:
# 1, intersection of r/t/n, 2, t/n/ 3. r/t, 4, r/n, 5 n only, 6 t only, 7 r only,
# output group 7
from generalUtils import *
import sys,getopt
import re

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hn:t:r:o:",["ifile=","ofile="])
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-n","--normal"):
        nf = arg
    elif opt in ("-t","--tumor"):
        tf = arg
    elif opt in ("-r","--relapse"):
        rf = arg
    elif opt in ("-o","--ofile"):
        output = arg
        outputCount = output + ".countIntersection"
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' +nf + "\t" + tf + "\t" + rf)
print('Output file:\t'+ output)


def getMut(filename):
    data = []
    cnt = 0 
    with open(filename) as f:
        line = f.readline()
        while line:
            if re.match("^#",line):
                pass
            else :
                tempChrom, tempPos, tempInfo = line.strip().split("\t",2)
                tempMut = Mutation(tempChrom,tempPos) 
                tempMut.info = tempInfo 
                # if cnt < 10:
                    # tempMut.disp()
                data.append(tempMut) 
                cnt = cnt + 1
            line = f.readline()
    print str(cnt) + "\t mutations processed" 
    return data

relapse = getMut(rf)
tumor = getMut(tf)
normal = getMut(nf)

allDict = {}

for mut in normal: 
    allDict[(mut.chr,mut.ps)] = 1

for mut in tumor :
    if allDict.get((mut.chr,mut.ps),0):
        allDict[(mut.chr,mut.ps)] = allDict[(mut.chr,mut.ps)] + 2
    else:
        allDict[(mut.chr,mut.ps)] = 2

for mut in relapse :
    if allDict.get((mut.chr,mut.ps),0):
        allDict[(mut.chr,mut.ps)] = allDict[(mut.chr,mut.ps)] + 4
    else:
        allDict[(mut.chr,mut.ps)] = 4

cnt1 = 0; cnt2 = 0; cnt3 = 0; 
cnt4 = 0; cnt5 = 0; cnt6 = 0; 
cnt7 = 0; 
outputH = open(output,'w') 
outputCountH = open(outputCount, 'w')

for k, v in allDict.items():
    if v == 1:
        cnt1 = cnt1 + 1
    elif v == 2:
        cnt2 = cnt2 + 1
    elif v == 3:
        cnt3 = cnt3 + 1
    elif v == 4:
        cnt4 = cnt4 + 1
        info = [mut.info for mut in relapse if (mut.chr, mut.ps) == k]
        outRec = "\t".join(map(str, list(k))  + info) #+"\t"+  k.info + "\n"
        outputH.write(outRec + "\n")
    elif v == 5:
        cnt5 = cnt5 + 1
    elif v == 6:
        cnt6 = cnt6 + 1
    elif v == 7:
        cnt7 = cnt7 + 1
    else:
        pass

outputCountH.write("normal\t" + "tumor\t" +  "relapse\t" +  "norm_tum\t" + \
"norm_relap\t" +  "tum_relap\t" +   "ntr\n") 
outputCountH.write( str(cnt1) + "\t" +  str(cnt2) + "\t" \
                   + str(cnt4) + "\t" +  str(cnt3) + "\t" \
                   + str(cnt5) + "\t" + str(cnt6) + "\t"\
                   + str(cnt7) + "\n")
