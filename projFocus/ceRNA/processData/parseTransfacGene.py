#!/usr/bin/python
#J.HE
'''Desp. : pasrse TRANSFAC gene.txt file, extract required information which will
 used to group mutations. Accesion number(AC); Gene name; DE(description);
 SY(symbol); CH(chrom); BS(binding site); BR(binding region)  '''

import sys,getopt
import re, os 
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
        outlog = output + ".log"
print('Script path'+ sys.argv[0])
print('Input file:' + input)
# print('Output file:'+ output)

def myreadlines(f, newline):
    buf = ""
    while True:
        while newline in buf:
            pos = buf.index(newline)
            yield buf[:pos]
            buf = buf[pos + len(newline):]
        chunk = f.read(4096)
        if not chunk:
            yield buf
            break
        buf += chunk

def getField(line,head):
    pattern = "^"+head
    if re.match(pattern,line):
        return line.strip().split()[1]
    else:
        pass


cnt = 0
cntBS = 0 
outfH = open(output,'w')
with open(input) as f:
    for record in myreadlines(f, "//"):
        # print record
        gene = {}
        for line in record.split("\n"):
            # print line
            # if re.match("^AC", line):
                # print line.strip().split()[1]
                # gene['AC'] = line.strip().split()[1]
            if re.match("^ID", line):
                temp = re.split("[\s-]*|\$",line.strip())
                species, gs = re.split("\$",temp[1])
                if species == "HS":
                    cnt = cnt + 1
                    outfH.write(gs+"\n")                    
                # if re.match("^OS", line) and re.find("Human",line.strip()):
                    # break
            # if re.match("^SD", line):
            #     gene['SD'] = line.strip().split()[1]
            #     print line.strip().split()[1]
            #     if re.match("^BS", line) and re.find("HS\$", line):
            #         # bsloc, bsAC, bsID, bsfact, bsfactorAC =\
            #         print re.split("\t|;",line.strip())
            #         if gene.get('BS',0):
            #             gene['BS'].append[()
            #         else:
            #             gene['BS'] = 1



            # print getField(line, "AC")
        # break
                
print "total genes \t " + str(cnt)    
outfH.close()
