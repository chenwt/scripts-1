import subprocess
import glob
import os
pid="PASFEW"
f1=pid + "-NoA.freq.temp"
f2=pid + "-TuA.freq.temp"
f3=pid + "-ReA.freq.temp"
print fn1, fn2, fn3
with open(f1) as input1:
  dataNo=zip(*(line.strip().split('\t') for line in input1))
with open(f2) as input2:
  dataTu=zip(*(line.strip().split('\t') for line in input2))
with open(f3) as input3:
  dataRe=zip(*(line.strip().split('\t') for line in input3))

