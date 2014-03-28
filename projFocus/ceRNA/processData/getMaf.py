#!/usr/bin/python
#J.HE
#Desp.: get som matrix for given gene list from a list of fils  
#input: 1.<VCF files directory> 2.<grandularity: snp/interval/gene > 3.< types of overlap with genes> 4. < gene list: first colun gene> 
#output: <file: mat>  
#update: with cutoff of up/down stream to 2kb 

import os,getopt
import re
import sys
import operator
from optparse import OptionParser
from subprocess import Popen,PIPE

usage = "python " + sys.argv[0] + " \n\
        -d <VCF file directory>  \n\
        -g <genelist, with exStart, exEnd, targeted> \n\
        -o <output .mat file> \n\
        -r <locType: gene/ tss> \n\
        -c <distCut: 2000 defalt>\n\
        -k <keyType: maf/mut/>"
example = "python " + sys.argv[0] + "\n\
        -d /Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/mafs/test\n\
        -o test.out \n\
        -l tss \n\
        -c 2000 \n\
        -k maf\n\
        -g /Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list.geneSingleStartEnd "

MSG = "MESSAGE:\t"
ERR = "ERROR:\t"
SUC = "SUCESS\t" 
DISTCUT = 2000

parser = OptionParser()
parser.add_option("-d", type="string", dest="vcfDir")
parser.add_option("-o", type="string", dest="output")
parser.add_option("-l", type="string", dest="locType")
parser.add_option("-c", type="int",    dest="DISTCUT")
parser.add_option("-g", type="string", dest="genelist")
parser.add_option("-k", type="string", dest="keyType")

options, args = parser.parse_args()
optionDict = vars(options)

vcfdir = optionDict['vcfDir'] 
output = optionDict['output'] 
genelist = optionDict['genelist'] 
keyType= optionDict['keyType'] 
locType = optionDict['locType'] 
DISTCUT = optionDict['DISTCUT']

class Mutation:
  ''' class store one record from VCF file, include: chrom, posStart, posEnd, gene1, gene2, DP, ADref, ADalt, '''
  def __init__(self, chrom, posStart, posEnd):
    self.chrom = chrom 
    self.posStart = long(posStart)
    self.posEnd  = long(posEnd)
    self.name = ":".join(map(str,[chrom,posStart,posEnd])) 
  def addVal(self, valueDict):
    self.vals = valueDict
  def dispMut(self,key):
    print "\t".join(map(str,[self.chrom,self.posStart, self.posEnd, self.vals[key]]) )
  def getMAF(self):
    key = 'DP4'
    valKey = self.vals.get(key,'NA').strip()
    tempDp4 = map(float,valKey.split(","))
    tempRef = tempDp4[1] + tempDp4[2] 
    tempAlt = tempDp4[-1] + tempDp4[-2] 
    if tempRef > 0:
        maf = tempAlt / (tempRef + tempAlt)
    elif tempRef == 0:
        maf = 1
    else:
        maf = 'NA'
    return maf        

class Gene:
  def __init__(self,name, chrom, posStart, posEnd, strand):
      self.name = name 
      self.chrom = chrom
      self.posStart = long(posStart)
      self.posEnd = long(posEnd)
      self.strand = strand
      self.tss = ":".join(map(str,[chrom,posStart,posEnd]))
  def hasMut(self, mutation):
      if self.name.lower() == mutation.vals.get('geneName','0').strip().lower():
          return 1
      else :
          return 0
  def nearMut(self, mutation):
      # distCut = DISTCUT
      distCut = 2000 
      if self.chrom == mutation.chrom:
          if self.posStart > mutation.posEnd + distCut or self.posEnd < mutation.posStart - distCut :
              return 0
          else :
              return 1
      else :
          return 0
  def dispGene(self):
      print "\t".join(map(str,\
      [self.name,self.chrom,self.posStart, self.posEnd, self.strand]))
 
def newMutation(vcfRecord): 
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  H_LS-BH-A0E0-11A-13D-A128-09
    recArray = vcfRecord.strip().split("\t") 
    som = Mutation(recArray[0], int(recArray[1]), int(recArray[1]) + len(recArray[4]) )
    infoDict = dict(item.split("=",1) for item in  recArray[7].split(";") if "=" in item)
    formatDict = dict(zip(recArray[-2].split(":"),recArray[-1].split(":"))) 
    allDict = dict(infoDict, **formatDict)
    som.addVal(allDict)
    return som
   
def sortedGeneList(geneDict):
    return(sorted(geneDict,key=operator.attrgetter('name')))

##----get all target genelist
geneDict = {}
with open(genelist) as f:
     for row in f.readlines() : 
         row = row.strip().split("\t")
         geneDict[Gene(row[0],row[1],row[2],row[3],row[4])] = ''
print "number of input genes: \t " + str(len(geneDict)) 

vcfs = [file for file in os.listdir(vcfdir) if re.findall(r"vcf$", file)]
cntFile = 1
keyType = keyType.lower()
locType = locType.lower()

os.chdir(vcfdir)

if keyType == 'mut' and locType == 'gene':  
    outputH = open(output,'w')
    for vcf in vcfs:
      print vcf
      geneDict = dict((k,0) for (k, v) in geneDict.items())
      with open(vcf, buffering = 10870) as f:
        for line in f.readlines():
	     if re.findall(r"^##", line):
		  continue
	     elif re.findall(r"^#", line):
		  tempBarcode = line.strip().split("\t")[-1]
	     else:
		  tempMut = newMutation(line)
		  for gene in geneDict.keys():
   	             if gene.hasMut(tempMut) :
   	                 geneDict[gene] = max(geneDict.get(gene,0),1)
   	                 continue ## one mutation only in one gene
   	             else :
   	                 continue
      geneList = sortedGeneList(geneDict) ##changed to list type              
      if cntFile == 1:
          gnames = [g.name for g in geneList]
          outputH.write("geneName\t" + "\t".join(sorted(gnames)) + "\n"  )
          valOut = tempBarcode
          for key in geneList:
               valOut = valOut + "\t" + str(geneDict[key])
          outputH.write(valOut + "\n" )
          outputH.flush()
      else :
          valOut = tempBarcode
          for key in geneList:
               valOut = valOut + "\t" + str(geneDict[key])
          outputH.write(valOut + "\n" )
          outputH.flush()
      cntFile = cntFile + 1   
    outputH.close() 	              
    print "num of vcf files:\t" + str(cntFile-1)  
    PYTHON = "~/tools/python/Python_current/python"
    SRCgetMutMat = "~/bin/trfile"
    output2 = output + ".mat" 
    pGenMat = Popen(  SRCgetMutMat +" " + output + " " + \
                    output2 , stdout=PIPE, shell=True)
    err = pGenMat.communicate()[1]
    if not err:  
        pRm = Popen("rm " + output, shell=True)
        pRm.communicate()
        print SUC 
    else:
        print ERR + "in Getting Matrix file.... "
	sys.exit()
elif keyType == 'maf' and locType == 'tss':
    outputH = open(output,'w')
    outputH.write("\t".join(["geneName","mutCode","MAF","barcode" ]) + "\n")
    cntTemp = 0
    cntVCF = 0
    chrList = map(str,range(1,23)) + ['X','Y','x','y'] 
    geneDictArray = [[] for _ in range(len(chrList) + 1) ]
    for gene in geneDict.keys():
        try :
            geneDictArray[chrList.index(gene.chrom)].append(gene)
        except:
            geneDictArray[-1].append(gene)
    for vcf in vcfs:
      print vcf
      cntVCF = cntVCF + 1
      with open(vcf, buffering = 100087000 ) as f:
        for line in f.readlines():
            if re.findall(r"^##", line):
	            continue
            elif re.findall(r"^#", line):
	            tempBarcode = line.strip().split("\t")[-1]
            else:
                cntFile = cntFile + 1   
                tempMut = newMutation(line)
                try:
                    for gene in geneDictArray[chrList.index(tempMut.chrom)]:
                        if gene.nearMut(tempMut):
                            cntTemp = cntTemp + 1
                            if cntTemp % 10000 == 0 :
                                print str(cntTemp) +" pairs tested.."
                            valOut = "\t".join(map(str,[gene.name, tempMut.name, \
                                                    tempMut.getMAF(),tempBarcode]) ) + "\n"
                            outputH.write(valOut)
                        else:
                            pass
                except:
                    pass
    print "num of vcf files:\t" + str(cntVCF)  
    print "num of mutations:\t" + str(cntFile)
    print "num of output gene-mutation:\t" + str(cntTemp)
    outputH.close() 	              
    sys.exit()
else:
    print MSG + locType + ":" + keyType +" combinations are Underdevelopment!"
    sys.exit()
print "##-----END---"

