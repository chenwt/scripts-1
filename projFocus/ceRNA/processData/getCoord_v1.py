#!/usr/bin/python



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

outputTssH = open(output+"_Tss.tsv",'w+')
output5UTRH = open(output+"_5pUTR.tsv",'w+')
output3UTRH = open(output+"_3pUTR.tsv",'w+')

with open(input) as f:
    line = f.readline()
    line = f.readline()
    while line:
        crtline = line.strip().split("\t")
        # refgName  = crtline[0]  
        chrom = crtline[1][3:]  
        strand  = crtline[2]  
        txs  = crtline[3]  
        txe  = crtline[4]  
        cdss  = crtline[5]  
        cdse  = crtline[6]  
        exons  = crtline[8]  
        geneSymbol = crtline[10]
        if strand == "+" :
            thrpUTR = [cdse, txe]
            fivpUTR = [txs, cdss]
            tss = txs
        elif strand == "-":
            thrpUTR = [txs, cdss]
            fivpUTR = [cdse, txe]
        else:
            fivpUTR = ['0', '0']
            thrpUTR = ['0', '0']
            tss = 0
        if geneSymbol:
            outputTssH.write("\t".join([geneSymbol, chrom, tss,\
                                        str(long(tss) + 1), strand]) + "\n")
            output5UTRH.write("\t".join([geneSymbol, chrom] + fivpUTR + [strand]) \
                                + "\n")
            output3UTRH.write("\t".join([geneSymbol, chrom] + thrpUTR + [strand]) \
                                + "\n")
        else:
            pass
        line = f.readline()
print "#-------END------"
