#!/usr/bin/python
#J.HE
#Desp.: Given gene name, gslist, normalized expression matrix of tumor and normal samples,
# cernet network, run rscript to identify candidate driver regulator which can
# explain the disregulation of input gint target genes
# output: target gene, candidate regulator, coefficients, r-squares, p-values  

import sys,getopt, re
import pickle
import generalUtils as gu
from collections import defaultdict
from subprocess import Popen, PIPE
import subprocess 

argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hg:a:l:t:n:c:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-g"):
        gene = arg
    elif opt in ("-a"):
        annof = arg
    elif opt in ("-l"):
        gslistf = arg
    elif opt in ("-t"):
        exptf = arg
    elif opt in ("-n"):
        expnf = arg
    elif opt in ("-c"):
        cernetf = arg
    elif opt in ("-o"):
        output = arg
print('Script path:\t"'+ sys.argv[0])
print('Input file:\t' + gene)
print('tumor file:\t' + exptf)
print('normal file:\t' + expnf)
print('output file:\t'+ output)

cernet = pickle.load(open(cernetf + ".pickle"))
print len(cernet._r)
regs = cernet.getReg(gene) ## targene is the first one 
geneAnno = pickle.load(open(annof + ".pickle")) 
print len(geneAnno)

regObjlist= [] 
for g in regs:
    if geneAnno.get(g, 0):
        chrom, tss, tse, strand = geneAnno[g]
        gObj = gu.GeneCoord(g, chrom, tss)
        gObj.tse = long(tse)
        gObj.strand = strand
        regObjlist.append(gObj)
    else:
        continue
print "number of all regulators:\t" + str(len(regs))
print "number of all regulators:\t" + str(len(regObjlist))

smps = []
## get samples
with open(gslistf) as f:
    line = f.readline()
    while line:
        if re.findall(gene, line.strip()):
            smps = re.split("\t|;",line.strip())[1:]
            break
        line = f.readline()
print "number of samples:\t" + str(len(smps))

## get expression matrix
def getExpDict(filename, annoDict):
    geneByChr = defaultdict(dict) 
    cnt = 1
    with open(filename) as f:
        line = f.readline()
        while line:
            if cnt == 1:
                allSamples = line.strip().split("\t")
                geneByChr['samples']={'all': allSamples}
            else:
                tmpGene,vals = line.strip().split("\t",1)
                try :
                    geneByChr[annoDict[tmpGene][0]].update(\
                            {tmpGene:vals.split("\t") })
                except IndexError:
                    pass 
                    # outputlogH.write(tmpGene + "\n")
            cnt = cnt + 1
            line = f.readline()
    return geneByChr    

expTumByChr  = getExpDict(exptf, geneAnno)  
outTumTemp  = output + "_reg_exp_tumor.temp" 
outputTumH = open(outTumTemp, 'w')

smpIndex = [ idx for (idx, val) in enumerate(expTumByChr['samples']['all']) \
            if val in smps  ]

print "sample number", len(smpIndex) 
outputTumH.write('gene' + "\t" + \
        "\t".join(map(expTumByChr['samples']['all'].__getitem__, smpIndex)) + "\n")

cntReg = 0 
for regObj in regObjlist:
    try:
        regTumExp = map(expTumByChr[regObj.chr][regObj.name].__getitem__, smpIndex)
        outputTumH.write(regObj.name + "\t" + "\t".join(regTumExp) + "\n")
        cntReg = cntReg + 1
    except KeyError:
        print regObj

### call Rscript

outputTumH.close()
cmd = "/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript \
    /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-2_regKeyRegulators.r\
    --vanilla --input "\
     +  outTumTemp + " --output " + output
print cmd 
# p = Popen(cmd, shell = True, stderr = PIPE, stdout = PIPE)
rtncode  = subprocess.call(cmd, shell = True)
# stdout, err = p.communicate()
if rtncode != 0 :
    print "Error in regression!"
    sys.exit()

print "#---END----"


