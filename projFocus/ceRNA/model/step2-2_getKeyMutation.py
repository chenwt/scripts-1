#!/usr/bin/python
#J.HE
#Desp.: given key regulators, group together all 3 prime utr/ promoter region
#run association to get driver mutations.

import sys,getopt
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' \n\
        -r <regulator file> \n \
        -c <copy number data> \n \
        -m <methylation data > \n\
        -s <somatic mutation data > \n\
        -g < gslist > \n \  
        -e < expression data> \n \
        -o <output file > '
example = 'python'  + sys.argv[0] + ' -i <input>  -o <output>'

try:
    opts,args = getopt.getopt(argv,"hi:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i"):
        input = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + input)
print('Output file:\t'+ output)


with(open(input)) as f:
    line = f.readline()
    while line:
        line = line.strip()
        line = f.readline()


