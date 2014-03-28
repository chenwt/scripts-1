#!/usr/bin/python
#J.HE

import sys,getopt
import re
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

print ('Script path'+ sys.argv[0])
print('Input file:' + input)
print('Output file:'+ output)

outputH = open(output, 'w')
with open(input) as f:
    line = f.readline()
    while line:
        if  re.match("^[gene|#]", line):
            line.strip().split("\t", 2)[2]

            outputH.write("
        line = f.readline()
