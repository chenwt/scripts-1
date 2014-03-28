#!/usr/bin/python
#J.HE
#Desp> given Gint name, return expression file of Gint, and it's Regulator
# seperator
##TODO: make this into a common code to extract givne rows and columns from a
# well formated file



###init and setup
gsListFile=

import sys,getopt,os
argv = sys.argv[1:]
input = ''
output = ''
usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'
try:
    opts,args = getopt.getopt(argv,"hg:o:",["ifile=","ofile="])
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-g","--ifile"):
        input = arg
    elif opt in ("-o","--ofile"):
        output = arg
print ('Script path:'+ sys.argv[0])
print('Input file:' + input)
print('Output file:'+ output)

