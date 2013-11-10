#!/usr/bin/env python
#input:<string:pid>
#output:<file1: comm mutation>
#	<file2: TuOnly mutation>
#	<file3: ReOnly mutations>

import sys

#for arg in sys.argv:
pid=sys.argv[1]
if pid != '':
  print pid
  inpT = pid + "-Tu.somatic.FINAL"
  inpR = pid + "-Re.somatic.FINAL"
  outT = pid + "_TuOnly.txt"
  out3 = pid + "_CommTuRe.txt"
  outR = pid + "_ReOnly.txt"
else:
  print "ERROR!"
  sys.exit(2)
chrPosArray = []
valArray = []
with open(inpT) as inpf1:
  for line in inpf1.readlines():
    [chrom,pos,val]=line.strip().split("\t",2)
    chrPosArray.append(chrom + "_" + pos)
    valArray.append(val)

outfT = open(outT, 'w')
outfR = open(outR, 'w')
outf3 = open(out3, 'w')

with open(inpR) as inpf2:
  for i,line in enumerate(inpf2):
    if i == 0:
      [chrom,pos,val] = line.strip().split("\t",2)
      tempKey = chrom + "_" + pos
      header = tempKey + "\t" + val
      outfR.write(header + "\n")
    else:
      [chrom,pos,val] = line.strip().split("\t",2)
      tempKey = chrom + "_" + pos
      if tempKey in chrPosArray:
        idx = chrPosArray.index(tempKey)
        outf3.write(chrPosArray[idx] + "\t" + valArray[idx] + \
            "\t|\t" + val +"\n")
        chrPosArray.pop(idx)
        valArray.pop(idx)
      else:
        outfR.write(chrom + "_" + pos + "\t" + val + "\n")
outfT.write(header + "\n")
for elemt in zip(chrPosArray, valArray) :
  outfT.write("\t".join(elemt) + "\n")
print "----END----"
outfT.close()
outfR.close()
outf3.close()

