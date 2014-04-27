import re
def parseKeyRegFile(file, pval_cut):
    with(open(file)) as f:
        line = f.readline()
        if re.findall("/", file):
            tgene = file.strip().split("/")[-1]
        else :
            tgene = file 
        tgene = tgene.split("_",1)[0]
        tsum = [tgene, 0.0, 0.0, 0, 0]
        regs = []
        for i, line in enumerate(f):
            if re.match(r"^#r2\t", line):
               tsum[1] = line.strip().split("\t")[1] 
            elif re.findall(r"^#r2.pval", line):
                _, pval = line.strip().split("\t")
                if float(pval) <= pval_cut:
                    tsum[2] = str(pval)
                else:
                    tsum = pval 
                    regs = ''
                    break
            elif re.match(r"^#totalReg", line):
                tsum[3] = line.strip().split("\t")[1]
            elif re.match(r"^#sigReg",line):
                tsum[4] = line.strip().split("\t")[1]
            elif re.match(r"^#target", line):
                if not line.strip().split("\t")[1] == tgene:
                    print "error at gene name matching"
            elif re.match(r"^#regulator", line):
                pass 
            else:
                reg, beta, pval = line.strip().split("\t")
                if float(beta) > 0.:
                    regs.append(reg)
    # print tsum, regs 
    return tsum, regs 


