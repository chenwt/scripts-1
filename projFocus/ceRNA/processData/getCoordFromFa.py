#!/usr/bin/python
##J.HE
#ipnut: <refseq fasta file, downloaded from UCSC>
#output: gene with location, 3pad,5pad, each sequence header per line 
import re
import sys,getopt
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

outputH = open(output,'w+')
cnt = 0 
with open(input,buffering=10000000) as f:
    line = f.readline()
    while line:
        if re.match(r'^>',line):
            pt = "\s+|="
            g, range, rv, p5, p5v, p3, p3v, s, sv, other= re.split(pt, line.strip(),9)
            rc,rs,rend = re.split(':|-',rv)
            outputH.write("\t".join([g.split("_")[2], rc[3:], rs, rend, p5v, p3v, sv]) + "\n") 
            outputH.flush()
            cnt = cnt + 1
            if cnt % 10000 == 0:
                print str(cnt) + "\t line writed..."
        line = f.readline()
outputH.close()
print str(cnt) + "\tline processed"
print "##----END------"

