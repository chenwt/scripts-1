#!/usr/bin/python
#J.HE
# This is to use greedy methods to identify a subset of ceRNA regulator * mutated samples for one target cancer gene. 
# This subset should optimized the correlation between sum-up ceRNA driver expression and target expression.
# steps involved:

# 'Column method' to optimized differential expression(K-S test)
# -1 prepareData, extract all sample, make sure the sample order is the sample
# -2 calculated fullMatrix test statistics, p-value corr_0,
# -3 set tk-1 = corr_prev, flip one mutation from 1 to 0(deactive it)
# -4 calculate corr_k, fisher transformed, 
# -5 if corr_k - corr_k-1 > tolenrance: deactivate mutation , otherwise filp it back to 1 (active mutation) 
# -6 repeat step4-5 until all mutations have been visited.
# 
# 'select optimal tolenrance'
# The idea: the more stringent tol produce less mutation.
# select the tolerance which produe maxi correlation, with significant p-value compare to full matrix and permut matrix
# 
# 'select optimal over all random initiation' 
# The idea: random initiation will proceduce slight different result, so try to approximate to global optmization by randomization.
# -1 Ran multiply random initization ( n = 100), select the one with optimized correlation, with qualify p-values
# -2 count the number of mutations being selected, output the final mutation

# keyRegSumfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/keyRegSumfile.cancergene.05052014"
# expfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
# gslistfile="//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list"
# mutfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero"
# output="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014"
# logDir="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/log"

import  sys,getopt
from    collections     import defaultdict
from    parseKeyRegFile import parseKeyRegFile
from    collections     import Counter, Sequence
from    parseGslistFile import parseGslistFile
import  numpy           as np
import  subprocess
import  os.path

usage = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
example = 'python ' + sys.argv[0] + \
        '-n <tf network> \
        -g <target gene list> \
        -m <mutfile>\
        -e <expfile> \
        -o <output> \
        -t <tolenrence type, fix/flex> \
        -r <random init iteration, 1000/> \
        -p <random permut actual matrix, 100/> \
        -s <selection type, max/all/<any number < r> > \
        -l <logdir> '

cntperm = 100
nrand = 1000

argv = sys.argv[1:]
try:
    opts,args = getopt.getopt(argv,"hn:g:m:e:t:s:r:p:l:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-n"):
        netfile = arg
    elif opt in ("-m"):
        mutfile = arg
    elif opt in ("-g"):
        tgfile = arg
    elif opt in ("-e"):
        expfile = arg
    elif opt in ("-t"):
        ttol = arg
    elif opt in ("-s"):
        tsel = arg
    elif opt in ("-r"):
        nrand = arg
    elif opt in ("-p"):
        cntperm = int(arg)
    elif opt in ("-l"):
        logDir = arg
    elif opt in ("-o"):
        output = arg

# print('Inputs:\n' + gslistfile + "\n" + keyRegSumfile + "\n" + \
#      mutfile + "\n" + logDir + "\n" + "ttol: " + ttol + "\n" +\
#      "tsel:" +  tsel + "\n" + "nrand:" +  nrand  )
# print('Output file:\t'+ output)

def formatSampleName(code19):
    if len(code19) >11:
        return code19[5:16].replace("-",".")
    else :
        return code19.replace("-", ".")

def loadTgs(genefile):
    resL = []
    with open(genefile) as f:
        line = f.readline()
        while line:
            resL.append(line.strip().split()[0])
            line = f.readline()
    return(resL)

def getMutExpSmp(expfile, mutfile):
    resL = []
    with open(expfile) as f:
        resL = f.readline().strip().split("\t")
        resL = map(formatSampleName, resL)
    with open(mutfile) as f:
        allMutSmp = f.readline().strip().split("\t")
        resL = [a for a in map(formatSampleName, allMutSmp) if a in resL]
    return resL


def loadExp(expfile, smpsL):
    resD = defaultdict(list)
    with open(expfile) as f:
        allExpSmp = map(formatSampleName, f.readline().strip().split("\t"))

        smpIndx = [id for (id, v) in enumerate(allExpSmp) if v in smpsL]
		## sort sample order according to smpsL
		smpIndx = [ id for id,v in sorted(zip(smpsL, smpIndx))]

        resD['gene'] = map(allExpSmp.__getitem__, smpIndx)

        line = f.readline()
        while line:
            crt_g, crt_e = line.strip().split("\t",1)
            temp = map(float,crt_e.split("\t"))
            resD[crt_g] = map(temp.__getitem__, smpIndx)

            line = f.readline()
    return resD

def loadMut(mutfile, smpsL):
    resD = defaultdict(list)
    with open(mutfile) as f:
        _, allMutSmp = f.readline().strip().split("\t",1)
        allMutSmp = map(formatSampleName, allMutSmp.split("\t"))
        smpIndx = [id for (id, v) in enumerate(allMutSmp) if v in smpsL]
		smpIndx = [ id for id,v in sorted(zip(smpsL, smpIndx))]
        resD['gene'] = map(allMutSmp.__getitem__, smpIndx)
        line = f.readline()
        while line:
            crt_g, crt_v = line.strip().split("\t",1)
            temp =  map(int,crt_v.split("\t"))
            resD[crt_g] = map(temp.__getitem__, smpIndx)
            line = f.readline()
    return resD

def format2Dlist(l2d):
    out = []
    for i in l2d:
        out.append("["+",".join(map(str,i)) + "]")
    return ";".join(out)

def loadTarReg(net, tgL):
	resD = defaultdict(list)
	with (open(netfile) ) as f:
		line = f.readline()
		while line:
			g1, g2 = line.split()[:2]
			if g2 in tgL:
				resD[g2].append(g1)
			line = f.readline()
	return resD

### loading data
tgeneList = loadTgs(tgfile)
mutExpSmpL = getMutExpSmp(expfile, mutfile)

expD = loadExp(expfile, mutExpSmpL, tgeneList)

mutD = loadMut(mutfile, mutExpSmpL, tgeneList)

tarRegD = loadTarReg(netfile, tgeneList)

outputH = open(output + ".lt3keyReg", 'w') 
outputErrH = open(output + ".errTargetGene", 'w')

cntqsub = 0 

for tgene in tarRegD.keys():

    allRegsL = tarRegD[tgene]

    regMutD = {k:v for (k,v) in mutD.items() if k in allRegsL}
    regExpD = {k:v for (k,v) in expD.items() if k in allRegsL}
    
    expMutRegL = set(regExpD.keys()).intersection(set(regMutD.keys()))
    outTempMut = output + "_" + tgene + "_regMut.temp"
    outTempExp = output + "_" + tgene + "_exp.temp"

    if len(regExpD.keys()) <= 3:
        outputH.write(tgene + "\t" + ";".join(regMutD.keys()) + "\n")
        continue

    outTempMutH = open(outTempMut,'w')
    outTempExpH = open(outTempExp,'w')
    
    outTempMutH.write('gene\t'+"\t".join( map(mutD['gene'].__getitem__, intMutSmpIdL)) + "\n")
    for k,v in regMutD.items():
        if k in expMutRegL: 
            outTempMutH.write(k + "\t" + "\t".join(map(str,v)) + "\n")

    outTempExpH.write('gene\t'+"\t".join(map(expD['gene'].__getitem__, intExpSmpIdL)) + "\n")

    ###debug----
    # print tgene,len(tarRegD.keys())
    # print expD[tgene]
    # continue
    # ###debug--
    
    ## target expession in the first row
    outTempExpH.write(tgene + "\t" + "\t".join(map(str,map(expD[tgene].__getitem__, intExpSmpIdL))) +"\n")
    for k,v in regExpD.items():
        if k in expMutRegL:
            outTempExpH.write(k + "\t" + "\t".join(map(str,v)) + "\n")  
    outTempMutH.close()
    outTempExpH.close()
    
    outTemp = output + "_" + tgene + ".tsv"

    RSCRIPT = "/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript"
    LOGDIR = logDir
    JOBNAME = tgene + "_optCor"
    RCODE = "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.r"
    TTOL = ttol #"flex"
    TSEL = tsel #"max" 
    NRAND = nrand
    OUTPUT = output + "_" + tgene
    cmd = "qsub -l mem=8g,time=40:: -S /" +  RSCRIPT + \
            " -e " + LOGDIR + \
            " -o " + LOGDIR + \
            " -N " + JOBNAME + \
            " -cwd " + \
            RCODE + " --vanilla " + \
            " --exp " + outTempExp +\
            " --mut " + outTempMut + \
            " --output " + OUTPUT + \
            " --ttol "  + TTOL + \
            " --tsel " + TSEL + \
            " --nrand " + NRAND 
    # print cmd
    if not os.path.isfile(outTemp) : 
        subprocess.Popen(cmd, shell = True)
    else:
        print outTemp + " already exist, skip it " 

    for iperm in range(cntperm):
        JOBNAME = tgene + "_permu" + str(iperm)
        RCODE = "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr_permuAll.r"
        TTOL = ttol #"flex"
        TSEL = tsel #"max" 
        NRAND = nrand
        OUTPUT = output + "_" + tgene
        cmd = "qsub -l mem=8g,time=40:: -S /" +  RSCRIPT + \
                " -e " + LOGDIR + \
                " -o " + LOGDIR + \
                " -N " + JOBNAME + \
                " -cwd " + \
                RCODE + " --vanilla " + \
                " --exp " + outTempExp +\
                " --mut " + outTempMut + \
                " --output " + OUTPUT + "_" + str(iperm) + \
                " --ttol "  + TTOL + \
                " --tsel " + TSEL + \
                " --nrand " + NRAND 
        # print cmd
        if not os.path.isfile(outTemp) : 
            subprocess.Popen(cmd, shell = True)
        else:
            print outTemp + " already exist, skip it " 
    cntqsub = cntqsub + 1 + cntperm


print "summitted job: ", cntqsub 
outputH.close()
outputErrH.close()
print "#[END]"

