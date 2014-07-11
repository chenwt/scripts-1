#!/usr/bin/python
#J.HE
'''
input: 1. file name pattern
ouptput: pval, integrated correlation of permuation
'''

import os, sys

argv = sys.argv[1:]

fptn = argv[0]
output = argv[1]

fnpermA = sorted( [f for f in os.listdir(".") if \
       re.match(r'optCorr\.result.*flex_max_1000\.permuAll', f)] )
fnactuA = [f for f in os.listdir(".") \
           if not re.match(r'optCorr\.result_*flex_max_1000\.permuAll', f)\
           and re.match(r'optCorr\.result_*flex_max_1000',f)]

def getTarName(fnactu):
    return fnactu.split("_")[1]

tarA = map(getTarName, fnactuA)


