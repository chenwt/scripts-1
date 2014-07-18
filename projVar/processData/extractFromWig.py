#!/usr/bin/python
#J.HE
'''
input: -m mutation file <bed> 
        -f feature file < wig/bed>
output: -o bedfile with mutation 

'''

import sys,getopt
from collections import defaultdict

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hf:m:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-m"):
        mutfile = arg
    elif opt in ("-f"):
        fefile = arg
    elif opt in ("-o"):
        output = arg

print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + input)
print('Output file:\t'+ output)

mutDict = defaultdict(list)
chromL = map(str, range(1,25) )

with(open(mutfile)) as f:
    line = f.readline()
    while line:
        line = line.strip().split()
        
        line = f.readline()
