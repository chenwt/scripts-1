#!/usr/bin/python
#J.HE

import sys,getopt
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
example = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hg:m:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-g"):
        gfile = arg
    elif opt in ("-m"):
        mutfile = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t'+ sys.argv[0])
print('Input file:\t' + gfile + "\n" + "\t" + mutfile) 
print('Output file:\t'+ output )


geneList = []
with(open(gfile)) as f:
    line = f.readline()
    while line:
        gene = line.strip()
        geneList.append(gene)
        line = f.readline()
print "%s of gene inputted" % len(geneList)

outputH = open(output, 'wt')
cnt = 0 
with(open(mutfile)) as f:
    line = f.readline()
    outputH.write(line) 
    line = f.readline()
    while line:
        gene, infos = line.strip().split("\t",1)
        if gene in geneList:
            outputH.write(line)
            cnt = cnt + 1
        line = f.readline()

print "%s of mutations output" % cnt
