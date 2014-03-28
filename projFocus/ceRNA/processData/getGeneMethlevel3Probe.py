#!/usr/bin/python
#J.HE


import os, re
import sys,getopt
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hg:p:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-g","--ifile"):
        gfile = arg
    elif opt in ("-p","--ifile"):
        pfile = arg
    elif opt in ("-o","--ofile"):
        output = arg
        outlog = output + ".log"
print('Script path:\t"'+ sys.argv[0])
print('Input file:\t' + input)

glist = []
outputH = open(output, 'w')
with open(gfile) as f:
    line = f.readline()
    while line:
        glist.append(line.strip())
        line = f.readline()
print "gene number\t" + str(len(glist))

cnt = 0 
with open(pfile) as f:
    line = f.readline()
    while line:
        if re.match("^cg",line):
            probe, val, geneName, chrom, ps = line.strip().split("\t")
            if re.findall(";", geneName):
                for g in re.split(";", geneName):
                    if g in glist:
                        outputH.write(probe +"\t" + geneName + "\t" + chrom +\
                        "\t" + ps +  "\n" ) 
                        cnt = cnt + 1
                        break
                    else:
                        continue
            else:
                if geneName in glist:
                    outputH.write(probe + "\t" + geneName + "\t" + chrom + \
                                  "\t" + ps + "\n" )
                    cnt = cnt + 1
                else:
                    pass
        line = f.readline()
print "output probe:\t" + str(cnt)  

outputH.close()
