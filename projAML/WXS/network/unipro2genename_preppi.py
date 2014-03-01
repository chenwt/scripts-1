#!/usr/bin/env python
# J.HE
# TODO: add desp. input, output comment
# 		make it command line useable!
import os
import getopt 
import sys 
import pprint

targetfile='temp_preppi_int_3col.txt';
mappingFile='temp_unipro2genenamePrePPI.txt';
outputFile='preppi_int_3col_genename.txt';
# targetfile='test1.txt';
# mappingFile='test2.txt';
# outputFile='test1_mapped.txt';
print targetfile 
print mappingFile
print outputFile 
# # get command line opts
# config = { 
#     "filei":"file1.txt", 
#     "filem":"file2.txt", 
#     "fileo":"file3.txt", 
     
# }
# opts, args = getopt.getopt(sys.argv[1:], 'hi:m:o:',  
#       [ 
#         'filei=',  
#         'filem=',  
#         'fileo='
#         ] 
#       )
# for option, value in opts: 
#     if  option in ["-h","--help"]: 
#         print """ 
#         Usage: unipro2genename_preppi.py -i test1.txt -m test2.txt -o test1_mapped.txt
#         """ 
#     elif option in ['-i','--filei']: 
#         config["filei"] = value
#     elif option in ['-m','--filem']: 
#         config["filem"] = value
#     elif option in ['-o','--fileo']: 
#         config["fileo"] = value
# # print '\t'.join(config)
# # print '\n'
# pp = pprint.PrettyPrinter(indent=1)
# pp.pprint(config)

# creat a dict
dmap = {}
# with open(config["filem"]) as f:
cnt = 1
with open(mappingFile) as f:
     for line in f:
             lline = line.split()
             if (len(lline) == 2):
                     (key, val) = lline
                     dmap[key] = val
                     cnt = cnt + 1
             # if(cnt == 200):
             #         break

f.close();
# replace target file using mapping

fout = open(outputFile,'w',1);
cnt = 1;
with open(targetfile) as f:
	for line in f:
		crtline = line.split()
		if(dmap.has_key(crtline[0]) & dmap.has_key(crtline[1])):
			crtline[0] = dmap[crtline[0]]
			crtline[1] = dmap[crtline[1]]
			fout.write('\t'.join(crtline))
			fout.write("\n")
			cnt = cnt + 1

fout.close();
print cnt;

