#!/usr/bin/python
#J.HE
#Desp.: given a gene sample list(target cancer genes), output the Gint freq,
# Gint related regulator number, gslist with more than 1 regulator


import re
import sys,getopt
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
        cernet = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t"'+ sys.argv[0])
print('Input file:\t' + input)
print('Input file 2 :\t' + cernet)
print('Output file:\t'+ output)

#-----------------funcS
def swap(xx):
    if xx[0] > xx[1]:
        return([xx[1], xx[0]])
    else:
        return([xx[0], xx[1]])

def sortCernet(net):
    net=map(swap, net)
    return net.sort()

def countReg(g, cernetSorted):
    cnt = 0 
    for [r1, r2] in cernetSorted:
        if r1 == g or r2 ==g :
            cnt = cnt + 1
        else:
            continue
    return(cnt)

#-----------------funcE

regulator = []
with open(cernet) as f:
    line = f.readline()
    while line:
        r1, r2 = line.strip().split("\t",2)[:2]
        regulator.append([r1,r2])
        line = f.readline()

sortCernet(regulator)

outputH1 = open(output + ".GintSmpCount", 'w')
outputH2 = open(output + ".GintRegCount", 'w')
outputH3 = open(input + ".hasReg.list", 'w')
outputH1.write("tGene\tgintSampleCount\n")
outputH2.write("tGene\tgintRegCount\n")
outputH3.write("tGene\tgintSamples\n")
with open(input) as f:
    line = f.readline()
    while line:
        g, s = re.split("\t", line.strip())
        smp  = re.split(";", s) 
        cntR = countReg(g, regulator)
        outputH1.write(g + "\t" + str(len(smp)) + "\n")
        outputH2.write(g + "\t" + str(cntR) + "\n")
        if cntR > 0:
            outputH3.write(line)
        line = f.readline()

outputH1.close()
outputH2.close()
outputH3.close()
print "#-----------------END---"
