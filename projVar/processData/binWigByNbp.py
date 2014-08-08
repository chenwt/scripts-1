#!/usr/bin/python
#J.HE
'''
To prepare genomic features,
Bin conservation scores by N, take the maximum score as the binned ones
'''

# inwig = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/genomicFeature/conservation/test/test.wig"
# n = 100
# flushCut = 1000

import sys,getopt

argv = sys.argv[1:]
usage = 'python ' + sys.argv[0] + '  -i <input wig file for one chrom> \
    -b <bin size> -c < flush out cut size, lines>'
example = 'python ' + sys.argv[0] + ' \
    -i /ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/genomicFeature/conservation/test/test.wig \
    -b 10 -c 20 '

try:
    opts,args = getopt.getopt(argv,"hi:b:c:")
except getopt.GetoptError as e: 
    print "ERROR", e
    print usage + "\n" + example 
    sys.exit(2)

for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i"):
        inwig = arg
    elif opt in ("-b"):
        n = int(arg)
    elif opt in ("-c"):
        flushCut = int(arg)

outfile = inwig + "_" + str(n)  + ".bed" 

def flushOut(outblock, outfh) :
    outfh.write("\n".join(outblock) + "\n")
    outfh.flush()


header = ''
outfh = open(outfile, 'aw') ## append to existing file

with(open(inwig)) as f: 
    for i in range(6):
        header = header + f.readline()
    line = f.readline()
    
#     outfh.write(header + "\t" + line +  "chrom\tstart\tend\tvalue\n")
    stepType, chrom, span = line.strip().split()

    span = int(span.replace("span=","")); chrom = chrom.replace("chrom=","")
    
    line = f.readline()
    outblock = [] * flushCut
    ps_prev, val_prev = line.strip().split()
    ps_prev = int(ps_prev); val_prev = float(val_prev)

    while line:
        ps_crt, val_crt = line.strip().split()
        ps_crt = int(ps_crt); val_crt = float(val_crt)
        if  ps_crt - ps_prev < n:
            #within span
            
            #print ps_prev, ps_crt
            val_prev = max(val_crt, val_prev)
            #print "ifs"
            
        else:
            # outside span
            outblock.append("\t".join(map(str,[chrom, ps_prev, ps_crt, val_prev])))
            
            if len(outblock) >= flushCut:
                flushOut(outblock, outfh)
                outblock = [] * flushCut
                #print "output"
                
            ps_prev = ps_crt; val_prev = val_crt 
            #print "reset"
    
        line = f.readline()
        
    flushOut(outblock, outfh)
    
outfh.close()   

print "finish"    
