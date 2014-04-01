#!/usr/bin/python
#J.HE
#Desp.: given the output file from step2-1_getKeyReg.r
#run association to get driver mutations.


import sys,getopt
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' \n\
        -g <target gene>\n \
        -d <dir for results from> \n \
        -m <mutation file> \n\
        -b <binding site> \n\
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


