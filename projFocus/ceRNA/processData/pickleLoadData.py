#!/usr/bin/python
#J.HE
#Desp.: save data as python binary file, serilizing python processing.

usage = 'python " + sys.argv[0] + " -i <input>  -o <output>'
example = 'python " + sys.argv[0] + " -i <input>  -o <output>'

import sys, os, getopt
import pickle 
input = sys.argv[1:][0]
print input
output = input + ".pickle"

inputH = open(input) 
data = pickle.load(inputH) 

print "#----END------"
