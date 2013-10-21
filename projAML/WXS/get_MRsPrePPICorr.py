#!/usr/bin/env python
# J.HE
# input: 1 MRs gene;
#		 a file has mutated gene names, one per line
#		 PrePPI database with genename as interactor.
# output: MRs, first order interacting proteins, mutated gene
# TODO:


# file names
mrf = 'MRs.sorted.genes'
# mrf = 'testMRs.genes'
mtgf = 'result_varFreq_v1_overlap.txt.genes'
preppi='/ifs/scratch/c2b2/ac_lab/jh3283/database/preppi/preppi_int_3col_genename.txt_600_90' 

# get preppi data 
dppi_a = {}
dppi_b = {}
with open(preppi) as f:
	for line in f:
		[p1, p2, val] = line.split()
		if( not dppi_a.has_key(p1)):
			dppi_a[p1] = [p2]
		else :
			dppi_a[p1] = dppi_a[p1] + [p2]
		if( not dppi_b.has_key(p2)):
			dppi_b[p2] = [p1]
		else:
			dppi_b[p2] = dppi_b[p2] + [p1]

f.close();



# get mutGene correlated protein
def intersect(a, b): #function to get intersection of two lists
	if ((len(a) > 0) & (len(b) > 0)):
		res = list(set(a) & set(b))
	else:
		res=[]
	return res

outfile = 'MR_mutgene_preppi.txt' 
outfile2 = 'MR_mutgene_preppi_ALL.txt' 

fout = open(outfile,'w') 
fout2 = open('MR_mutGene_preppi_stat.txt','w')
fout3 = open(outfile2,'w') 

with open(mtgf) as fin:
	for line in fin: # get MRs first order correlated protein
		[mtg] = line.split()
		mtg_p = []
		if (dppi_a.has_key(mtg)):
			mtg_p = dppi_a[mtg]
		if (dppi_b.has_key(mtg)):
			mtg_p = mtg_p + dppi_b[mtg]
		if (len(mtg_p) > 0): 
			# print(str(len(mr_p)))
			with open(mrf) as f:
				flag =1 
				for mline in f:
					mr_p =[]
					[mr] = mline.split()
					if (dppi_a.has_key(mr)):
						mr_p = dppi_a[mr]
					if (dppi_b.has_key(mr)):
						mr_p = mr_p + dppi_b[mr]
					both = intersect(mr_p,mtg_p)
					fout2.write(mr + '\t' + str(len(mr_p)) + '\t' + mtg +
					 '\t' + str(len(mtg_p)) + '\tboth\t' +str(len(both))+ '\n')
					fout3.write(mr + ':' + ','.join(mr_p)  + ':' +  mtg + ':' + ','.join(mtg_p)  + ':' + ','.join(both) )
					if(len(mr_p) > 0 ):
						fout.write(mr + '\t'+ mtg + '\t' + ','.join(both) + '\n' )

fout.close()
fout2.close()
fout3.close()

				# 
						# fout.write(mr + '\t'+ mtg + '\t' + ','.join(intersect(mr_p,mtg_p)) + '\n' )
					# for key in dppi.keys():
					# 	if(mtg in key & float(dppi[key] > 5400)):
							# mtg_crtpp = mtg_crtpp + key
						# fout.write(mtg + '\t' + ','.join(intersect(mr_p,mtg_crtpp)) + '\t' + mtg )

 
				
