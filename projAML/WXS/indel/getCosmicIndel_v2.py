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
print('cosmic file:\t' + cosmic)
print('Output file:\t'+ output)

cnt = 0
with(open(cosmic)) as f:
    line = f.readline()
    cosmicDict = defaultdict(list)
    while line:
        if not re.match("^Mutation", line):
            tchrom, tps, tpe, tgene, ttype = line.strip().split("\t",4)
            tmut = gu.Mutation(tchrom, tps)
            tmut.pe = int(tpe)
            tmut.gene = tgene 
            try :
                tidx = re.search(r'ins|del|A|T|G|C|>',ttype)
                if tidx:
                    tidx = tidx.start()
                    tmut.type = ttype[tidx:]
            except AttributeError:
                continue 
            cosmicDict[tchrom].append(tmut) 
            cnt = cnt + 1
        line = f.readline()
print "cosmic mutation number\t", cnt 
outputH = open(output,'w')

cnt = 0
cntd = 0 
cntk = 0
with(open(input)) as f:
    line = f.readline()
    while line:
        if not re.match("^#",line):
            tchr, tps, _, tseqref, tseqalt, _, tfilter, tinfo, _,_ =\
                    line.strip().split("\t")
            if len(tseqref) > 15 or len(tseqalt) > 15:
                line = f.readline()
                cntd = cntd + 1
                continue
            if not tfilter == 'PASS':
                line = f.readline()
                cntd = cntd + 1
                continue
            else:
                if len(tseqref) > len(tseqalt):
                    type = "del" + tseqref.replace(tseqalt,'',len(tseqalt))
                elif len(tseqref) == len(tseqalt):
                    type = tseqref + ">" + tseqalt 
                else:
                    type = "ins" + tseqalt.replace(tseqref,'',len(tseqref))
                crtmut = gu.Mutation(tchr, tps)
                crtmut.pe = crtmut.ps + len(tseqref) - 1
                crtmut.type = type
                if  [x for x in cosmicDict[tchr] if x.ps == crtmut.ps and
                     x.type == crtmut.type]:
                    outputH.write(line)
                    cntk = cntk + 1
                else:
                    pass
        elif re.match("^#CHROM",line):
            outputH.write(line)
        else:
            pass
        cnt = cnt + 1
        line = f.readline()
print  "total mutation\t" + str(cnt) 
print "filtered mutation \t" + str(cntd) 
print  "cosmic mutation\t" + str(cntk) 

outputH.close()
