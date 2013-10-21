##----------------------------
#need to put this script to result folder, and also getReports in results folder
# create and open files
import subprocess
import glob
import os
subprocess.call(['./getReports.sh']) 
flist = sorted(glob.glob('*temp'))
fiter = iter(flist)
for fn1, fn2, fn3 in zip(fiter,fiter,fiter):
	print fn1, fn2, fn3
	#fn1 Re_N, fn2 Re_T, fn3 T_N
	with open(fn1) as f1:
	    lines1 = f1.read().splitlines()
	with open(fn2) as f2:
	    lines2 = f2.read().splitlines()
	with open(fn3) as f3:
	    lines3 = f3.read().splitlines()
##----------------------------
#list operations
	com13 = list(set(lines1) & set(lines3))
	union13 = list(set(lines1) | set(lines3))
#	diff13 = list(set(union13) - set(com13))
	diff13 = list(set(lines1) - set(com13))
	bothfind = list(set(lines2) & set(diff13))
	#fout=open(fn1[0:6] + '_stats.txt', 'w+')
	#fout2=open(fn1[0:6] + '_consistGenes.txt','w+')
	fout3=open(fn1[0:6] + '_overlapReN_TN.txt', 'w+')
##----------------------------
#output statistics
	#print >>fout, "PID\t", fn1[0:6]
	#print >>fout, "T_N\t", len(lines3)
	#print >>fout, "Re_N\t", len(lines1)
	#print >>fout, "Re_T\t", len(lines2)
	#print >>fout, "TN_ReN\t", len(com13)
	#print >>fout, "ReN_not_TN\t", len(diff13)
	#print >>fout, "consist\t",len(bothfind)
	#print >>fout2, '\n'.join(bothfind)
	print >>fout3, '\n'.join(com13)
	f1.close()
	f2.close()
	f3.close()
	#fout.close()
	#fout2.close()
	fout2.close()

##----------------------------
#delete files
for filename in flist:
	os.unlink(os.getcwd() + '/' + filename)

