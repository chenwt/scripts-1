#!/usr/bin/python
#J.HE
#input: 1.<VCF files directory> 2.<grandularity: snp/interval/gene > 3.< types of overlap with genes> 4. < gene list: first colun gene> 
#output: <file: mat>  

usage = 'python " + sys.argv[0] + " -d <VCF file directory>  -o <output .mat file> \
        -r <locType: gene; tss>  -g < genelist, with exStart, exEnd, targeted> -k <keyType: maf/DP4/AD/som>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'

import os,getopt
import re
import sys
import operator
#from optparse import OptionParser
#parser = OptionParser()
#parser.add_option("-d", type="string", dest="vcfDir")
#parser.add_option("-o", type="string", dest="output")
#parser.add_option("-l", type="string", dest="locType")
#parser.add_option("-g", type="string", dest="genelist")
#parser.add_option("-k", type="string", dest="keyType")

#options, args = parser.parse_args()
#optionDict = vars(options)
#print optionDict
#
#vcfdir = optionDict['vcfDir'] 
#output = optionDict['output'] 
#genelist = optionDict['genelist'] 
#keyType= optionDict['keyType'] 
#locType = optionDict['locType'] 

vcfdir = '/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/mafs/test'
output = 'test.out'
genelist = '/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list.geneSingleStartEnd'
keyType= 'snp' 
locType = 'gene'
 
print('Script path', sys.argv[0])
print('Input file:', vcfdir)
print('Output file:', output) 

class Mutation:
  ''' class store one record from VCF file, include: chrom, posStart, posEnd, gene1, gene2, DP, ADref, ADalt, '''
  def __init__(self, chrom, posStart, posEnd):
    self.chrom = chrom 
    self.posStart = posStart
    self.posEnd  = posEnd
  def addVal(self, valueDict):
    self.vals = valueDict
  def dispMut(self,key):
    print "\t".join(map(str,[self.chrom,self.posStart, self.posEnd, self.vals[key]]) )

class Gene:
  def __init__(self,name, chrom, posStart, posEnd, strand):
      self.name = name 
      self.chrom = chrom
      self.posStart = posStart
      self.posEnd = posEnd
      self.strand = strand
  def hasMut(self, mutation):
      if self.name.lower() == mutation.vals.get('geneName','0').strip().lower():
          return 1
      else :
          return 0
  def nearMut(self, mutation):
      if self.chrom == mutation.chrom  :
          pass ##<====UNDER WORK TODO
      else :
          return 0
  def dispGene(self):
      print "\t".join(map(str,\
      [self.name,self.chrom,self.posStart, self.posEnd, self.strand]))
#def getMutVal(mutation,key):
#  mutation.vals[key]  
#  return [som.chrom, som.posStart, som.posEnd,som.vals[key]] 

 
def newMutation(vcfRecord): 
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  H_LS-BH-A0E0-11A-13D-A128-09
  recArray = vcfRecord.strip().split("\t") 
  som = Mutation(recArray[0], int(recArray[1]), int(recArray[1]) + len(recArray[4]) )
  infoDict = dict(item.split("=",1) for item in  recArray[7].split(";") if "=" in item)
  formatDict = dict(zip(recArray[-2].split(":"),recArray[-1].split(":"))) 
  allDict = dict(infoDict, **formatDict)
  #print allDict
  som.addVal(allDict)
  return som
   
def sortedGeneList(geneDict):
    #geneDictSorted = sorted(geneDict, key=operator.attrgetter('name'))
    #return(geneDictSorted)
    return(sorted(geneDict,key=operator.attrgetter('name')))
    #return(sorted(geneDict,key=[k.name for k in geneDict.keys()]))

##----get all target genelist
geneDict = {}
with open(genelist) as f:
     for row in f.readlines() : 
         row = row.strip().split("\t")
         geneDict[Gene(row[0],row[1],row[2],row[3],row[4])] = ''
print "number of input genes: \t " + str(len(geneDict)) 

vcfs = [file for file in os.listdir(vcfdir) if re.findall(r"vcf$", file)]
#print "\t".join(vcfs) if len(vcfs) > 1 else vcfs
cntFile = 1
outputH = open(output,'w')
for vcf in vcfs:
  #geneDict.values() = ''
  geneDict = dict((k,0) for (k, v) in geneDict.items())
  with open(vcf) as f:
    for line in f.readlines():
      if re.findall(r"^##", line):
	  continue
      elif re.findall(r"^#", line):
	  tempBarcode = line.strip().split("\t")[-1]
      else:
	  tempMut = newMutation(line)
	  keyType = keyType.lower()
	  locType = locType.lower()
	  if keyType == 'snp' and locType == 'gene':
	      for gene in geneDict.keys():
	          if gene.hasMut(tempMut) :
	              geneDict[gene] = max(geneDict.get(gene,0),1)
	              continue ## one mutation only in one gene
	          else :
	              continue
	  elif keyType == 'maf' and locType == 'tss':
	      for gene in geneDict.keys():
	          if gene.nearMut(tempMut):
	              geneDict[gene].append(tempMut)
	          else :
	              pass
  geneList = sortedGeneList(geneDict) ##changed to list type              
  if cntFile == 1:
      gnames = [g.name for g in geneList]
      outputH.write("geneName\t" + "\t".join(sorted(gnames)) + "\n"  )
      valOut = tempBarcode
      for key in geneList:
          valOut = valOut + "\t" + str(geneDict[key])
      outputH.write(valOut + "\n" )
  else :
      valOut = tempBarcode
      for key in geneList:
          valOut = valOut + "\t" + str(geneDict[key])
      outputH.write(valOut + "\n" )
  cntFile = cntFile + 1   

print "num of vcf files:\t" + str(cntFile-1)  

outputH.close() 	              
	         # key = 'DP4'
          #        valKey = tempMut.vals.get(key,'NA').strip()
          #    #resKey = 1 if valKey in genelistDict.keys()
	         # tempMut.dispMut(key)
	         # tempDp4 = map(float,valKey.split(","))
	         # tempRef = tempDp4[1] + tempDp4[2] 
	         # tempAlt = tempDp4[-1] + tempDp4[-2] 
	         # if tempDp4[1] + tempDp4[2] > 0:
          #           valRes = tempAlt / (tempRef + tempAlt)
          #           print valRes
          #        else:
          #           valRes = 'NA'
          #           continue
          #elif keyType == 'snp':
          #    continue 
          #   

