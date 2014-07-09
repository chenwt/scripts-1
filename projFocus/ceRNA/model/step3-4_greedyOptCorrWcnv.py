#!/usr/bin/python
#J.HE
# This is to use greedy methods to identify a subset of ceRNA regulator * mutated samples for one target cancer gene. 
# This subset should optimized the correlation between sum-up ceRNA driver expression and target expression.
# steps involved:

# 'Cell method'
# -1 prepareData
# -2 calculated fullMatrix correlation corr_0, fisher transformed 
# -3 set corr_k-1 = corr_prev, flip one mutation from 1 to 0(deactive it)
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
from    pandas          import DataFrame
import  pandas          as pd

usage = 'python ' + sys.argv[0] + ' -i <input>  -o <output>'
example = 'python ' + sys.argv[0] + \
        '-i <gslistfile> \
        -k <keyreg sum file> \
        -m < mutfile>\
        -e <expfile> \
        -o <output> \
        -t <tolenrence type, fix/flex> \
        -r <randome init iteration, 1000/> \
        -s <selection type, max/all/<any number < r> > \
        -l <logdir> '


argv = sys.argv[1:]
try:
    opts,args = getopt.getopt(argv,"hi:k:m:c:e:t:s:r:l:o:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-i"):
        gslistfile = arg
    elif opt in ("-m"):
        mutfile = arg
    elif opt in ("-c"):
        cnvfile = arg
    elif opt in ("-k"):
        keyRegSumfile = arg
    elif opt in ("-e"):
        expfile = arg
    elif opt in ("-t"):
        ttol = arg
    elif opt in ("-s"):
        tsel = arg
    elif opt in ("-r"):
        nrand = arg
    elif opt in ("-l"):
        logDir = arg
    elif opt in ("-o"):
        output = arg
print('Inputs:\n' + gslistfile + "\n" + keyRegSumfile + "\n" + \
     mutfile + "\n" + logDir + "\n" + "ttol: " + ttol + "\n" +\
     "tsel:" +  tsel + "\n" + "nrand:" +  nrand  )
print('Output file:\t'+ output)

def formatSampleName(code19):
    if len(code19) >11:
        return code19[5:16].replace("-",".")
    else :
        return code19.replace("-", ".")

def loadTarReg(keyRegSumfile):
    resD = defaultdict(list)
    with open(keyRegSumfile) as f:
        line = f.readline()
        while line:
            crt_t, ctr_r = line.strip().split("\t")
            resD[crt_t] = ctr_r.split(";")
            line = f.readline()
    return(resD)

def getMutExpSmp(expfile, mutfile):
    resL = []
    with open(expfile) as f:
        resL = f.readline().strip().split("\t")
        resL = map(formatSampleName, resL)
    with open(mutfile) as f:
        allMutSmp = f.readline().strip().split("\t")
        resL = [a for a in map(formatSampleName, allMutSmp) if a in resL]
    return resL

def loadTarIntSmp(gslistfile):
    resD = defaultdict(list)
    with open(gslistfile) as f:
        line = f.readline()
        while line:
            crt_t, crt_s = line.strip().split("\t")
            resD[crt_t] = crt_s.split(";")
            line = f.readline()
    return(resD)

def loadExp(expfile, smpsL):
    resD = defaultdict(list)
    with open(expfile) as f:
        allExpSmp = map(formatSampleName, f.readline().strip().split("\t"))
        smpIndx = [id for (id, v) in enumerate(allExpSmp) if v in smpsL]
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
        resD['gene'] = map(allMutSmp.__getitem__, smpIndx)
        line = f.readline()
        while line:
            crt_g, crt_v = line.strip().split("\t",1)
            temp =  map(int,crt_v.split("\t"))
            resD[crt_g] = map(temp.__getitem__, smpIndx)
            line = f.readline()
    return resD

def __test__():
    from collections import defaultdict
    seqSet  = {'G1':['S2','S4','S6'],
                'G2':['S1','S3'],
                'G3':['S1'],
                'G4':['S1'],
                'G5':['S5'],
                'G6':['S3']}
    seq     = ['S1', 'S2', 'S3', 'S4', 'S5','S6']
    weightSet = {'G1':[1.0,0.5,1.5],
                 'G2':[2.0, 2.5],
                 'G3':[2.3],
                 'G4':[1.2],
                 'G5':[2.5],
                 'G6':[3.0]}
    setObjL = []
    for sk, ss in seqSet.items():
        setObjL.append(MutSet(sk,ss,weightSet[sk]))
    geneL, smpL, costL = wgsc(setObjL, seq, wtype = "mean")
    geneL, smpL, costL = wgsc(setObjL, seq, wtype = "total")
    geneL, smpL, costL =  wgsc(setObjL, seq, wtype = "max")

def format2Dlist(l2d):
    out = []
    for i in l2d:
        out.append("["+",".join(map(str,i)) + "]")
    return ";".join(out)

def loadCNV(cnvfile, smpsL):
    '''
    load cnv matrix; and convert to binary 1/0
    '''
    resD = defaultdict(list)
    with open(cnvfile) as f:
        _, allMutSmp = f.readline().strip().split("\t",1)
        allMutSmp = map(formatSampleName, allMutSmp.split("\t"))
        smpIndx = [id for (id, v) in enumerate(allMutSmp) if v in smpsL]
        resD['gene'] = map(allMutSmp.__getitem__, smpIndx)
        line = f.readline()
        while line:
            crt_g, crt_v = line.strip().split("\t",1)                                  
            temp =  map(lambda x: int(float(x) != 0), crt_v.split("\t"))
            resD[crt_g] = map(temp.__getitem__, smpIndx)
            line = f.readline()
    return resD

def jaccardCoeff(lista, listb):
    cnti = 0.0; cntj = 0.0; cntmatch = 0.0 
    if len(lista) == 0 and len(listb) == 0:
        return 0
    else:
        for i,j in zip(lista, listb):
            if i !=0:
                cnti = cnti + 1
                if j != 0:
                    cntj = cntj + 1
                    if i == j:
                        cntmatch = cntmatch + 1
    if (cnti + cntj) > 0 :    
        return cntmatch / (cnti + cntj)        
    else :
        return 0

def mergeMutCnv(mutD, cnvD):
    '''
        sum up mutation matrix with cnv matrix, binary sum up, 
        1-0, 0-1, 1-1 ==> 1; 0-0 ==> 0
        sample order follo
    '''
    mutSmp = mutD['gene']; cnvSmp = cnvD['gene']
    sample = list(set(mutSmp).intersection(set(cnvSmp)))
    del mutD['gene'];del cnvD['gene']

    
    mutGene = mutD.keys(); cnvGene = cnvD.keys()
    gene = list(set(mutGene).intersection(set(cnvGene)))
    mutOnlyGene = list(set(mutGene) - set(gene))
    cnvOnlyGene = list(set(cnvGene) - set(gene))
    
    mutDF = pd.DataFrame(mutD, index=mutSmp).T
    cnvDF = pd.DataFrame(cnvD, index=cnvSmp).T
    
    #should concatenating two list union of all genes
    vartDF = mutDF.loc[gene,sample] + cnvDF.loc[gene, sample]  
    vartDF = pd.concat([vartDF, mutDF.loc[mutOnlyGene,sample], \
                     cnvDF.loc[cnvOnlyGene, sample]])
    
    
#     print mutDF.loc[gene[100:105],sample[20:25]]
#     print cnvDF.loc[gene[100:105],sample[20:25]]
#     print vartDF.loc[gene[100:105],sample[20:25]]
    
    
    mutD['gene'] = mutSmp; cnvD['gene'] = cnvSmp
    
    vartDF = vartDF.where(vartDF == 0, 1)
    vartD = vartDF.T.to_dict(outtype='list')
    vartD['gene'] = sample
    return vartD





tarRegD = loadTarReg(keyRegSumfile)
tarIntSmpD = loadTarIntSmp(gslistfile)
mutExpsmpL = getMutExpSmp(expfile, mutfile)

cnvD = loadCNV(cnvfile, mutExpsmpL)
mutD = loadMut(mutfile, mutExpsmpL)
expD = loadExp(expfile, mutExpsmpL)
vartD = mergeMutCnv(mutD, cnvD)


outputH = open(output + ".lt3keyReg", 'w') 
outputErrH = open(output + ".errTargetGene", 'w')
cntqsub = 0 
for tgene in tarRegD.keys(): 
    tIntSmp = tarIntSmpD[tgene] 
    allRegsL = tarRegD[tgene]
    intMutSmpIdL = [id for (id, s) in enumerate(vartD['gene']) if s in tIntSmp]
    intExpSmpIdL = [id for (id, s) in enumerate(expD['gene']) if s in tIntSmp]
    
    print tgene 
    regMutD = {k:map(v.__getitem__, intMutSmpIdL) for (k,v) in vartD.items() if k in allRegsL}
    regExpD = {k:map(v.__getitem__, intExpSmpIdL) for (k,v) in expD.items() if k in allRegsL}
    
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
    # cmd =  RSCRIPT + " " + RCODE + " --vanilla " + \
    #         " --exp " + outTempExp +\
    #         " --mut " + outTempMut + \
    #         " --output " + OUTPUT + \
    #         " --ttol "  + TTOL + \
    #         " --tsel " + TSEL + \
    #         " --nrand " + NRAND 
   
    if not os.path.isfile(outTemp) : 
        subprocess.Popen(cmd, shell = True)
    else:
        print outTemp + " already exist, skip it " 
    cntqsub = cntqsub + 1

print "summitted job: ", cntqsub 
outputH.close()
outputErrH.close()
print "#[END]"

