#!/usr/bin/python
#J.HE
'''
change the e+007 notation to numbers
'''

import sys, os

input = sys.argv[1:][0]
output = input + ".fixed"

def sci2num(x):
    fc1, fc2 = x.replace("e","").split("+")
    return( int(float(fc1) * (10 ** int(fc2)) ) )

outf = open(output, 'w')

with(open(input)) as f:
    line = f.readline()
    while line:
        chrom, ps, pe, info = line.strip().split("\t")
        # print ps
        if ps[-5:-3]=="e+" :
            # print sci2num(ps)
            outf.write(chrom + "\t" + str(sci2num(ps)) + "\t"\
                       + str(sci2num(pe)) + "\t" + info + "\n") 
        else:
            outf.write(line)
        line = f.readline()

outf.close()


print '[----END----]'
