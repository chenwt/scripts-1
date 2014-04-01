#!/usr/bin/python
#J.HE

import sys,getopt
import generalUtils as gu
import re
from collections import defaultdict

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hi:c:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i"):
        input = arg
    elif opt in ("-c"):
        cosmic = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t"'+ sys.argv[0])
print('Input file:\t' + input)
print('Output file:\t'+ output)

cnt = 0
with(open(cosmic)) as f:
    line = f.readline()
    cosmicDict = defaultdict(list)
    while line:
        if re.findall(r'ind|del',line):
            tchrom, tps, tpe, tgene, ttype = line.strip().split("\t",4)
            tmut = gu.Mutation(tchrom, tps)
            tmut.pe = int(tpe)
            tmut.gene = tgene 
            tidx = re.search(r'ins|del',ttype).start()  
            if tidx >=0:
                tmut.type = ttype[tidx:]
            cosmicDict[tchrom].append(tmut) 
            cnt = cnt + 1
        line = f.readline()

print cnt 
outputH = open(output,'w')
cnt = 0
cntd = 0 
cntk = 0
with(open(input)) as f:
    line = f.readline()
    while line:
        if not re.match("^#",line):
            tchr, tps, _, tseqref, tseqalt, _, _, tinfo, _,_ =\
                    line.strip().split("\t")
            if len(tseqref) > 15 or len(tseqalt) > 15:
                line = f.readline()
                cntd = cntd + 1
                continue
            else:
                crtmut = gu.Mutation(tchr, tps)
                crtmut.pe = crtmut.ps + len(tseqref) - 1
                if  [x for x in cosmicDict[tchr] if x.ps == crtmut.ps]:
                    outputH.write(line)
                    cntk = cntk + 1
                else:
                    pass
        else:
            outputH.write(line)
        cnt = cnt + 1
        line = f.readline()
print input + "total mutation\t" + str(cnt) 
print input + "long mutation \t" + str(cntd) 
print input + "cosmic mutation\t" + str(cntk) 

outputH.close()
