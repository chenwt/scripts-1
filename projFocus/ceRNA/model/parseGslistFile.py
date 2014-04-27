import re
def parseGslistFile(gene, file):
    # smplist = [] 
    smps = ""
    with(open(file)) as f:
        line = f.readline()
        while line:
            ptn = gene + "\t"
            if re.match(ptn, line):
                smps  = line.strip().split("\t",1)[1].split(";")
                break 
            line = f.readline()
    return smps
