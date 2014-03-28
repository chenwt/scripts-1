#!/usr/bin/python
#J.HE
#Desp.: combine promoter, utr3 and utr5 mutations for each gene 

import re
import sys,getopt
import generalUtils as gu 

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hp:t:f:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-p","--ifile"):
        inputp = arg
    elif opt in ("-t","--ifile"):
        input3 = arg
    elif opt in ("-f","--ifile"):
        input5 = arg
    elif opt in ("-o","--ofile"):
        output = arg
print('Script path:\t"'+ sys.argv[0])
print('Input file:\t' + input)

outputH = open(output, 'w')
hflag = 0 
for file in [inputp, input3, input5]:
    print file
    with open(file) as f:
        tgdlist = [{} for _ in range(25)]  
        line = f.readline()
        while line:
            tgene, tchr, tps, tpe, mgene, mps, mpe, vals =\
                        line.strip().split("\t",7)
            if not re.findall('gene', line) :
                ntchr = gu.chr2Num(tchr)
                ttgdlist = tgdlist[ntchr]
                if ttgdlist.get(tgene,0) :
                    tgdlist[ntchr][tgene] = map(lambda x:x[0] + x[1],\
                            zip(ttgdlist.get(tgene), map(int, vals.split("\t")))) 
                else:
                    tgdlist[ntchr][tgene] = map(int,vals.split("\t"))
                    hflag = 1
            else:                 
                if hflag == 0:
                    outputH.write("targetGene" + "\t" + vals + "\n") 
                else:
                    pass 
            line = f.readline()

tgdlist = filter(lambda x:len(x)>0,tgdlist)
for i in range(len(tgdlist)):
    for key,val in tgdlist[i].items() :
        outputH.write(key + "\t" + "\t".join(map(str, val)) + "\n" ) 
outputH.close()

