#!/usr/bin/python
#input: <file:TCGA barcode each line>
#output: <file: TCGA barcode \tabe sample_code>
#usage: python do_designTCGAsample.py barcode.list

inf="barcode.list"
outf="barcode.list.sample_code"
fsamplecode = "/ifs/scratch/c2b2/ac_lab/jh3283/projAML/sampleInfo/sampleType.txt"
df = open(fsamplecode)
scDict = {} 
for i,line in enumerate(df):
	if i > 0:
		line = line.strip()
		[key,des,value] = line.split("\t")
		scDict[key] = value
	else :
		next
# print len(scDict.items())
# for k,v in scDict.items():
# 	print k

f = open(inf)
fout = open(outf,"w")
for i, line in enumerate(f):
	temp = line.strip()
	tempsc = temp.split("-")[3][:2]
	# print temp
	# print tempsc
	fout.write(temp + "\t" + scDict[tempsc] + "\t" + "pair-end" + "\n")

f.close()
fout.close()



