#!/usr/bin/python
#J.HE
#Desp.: save data as python binary file, serilizing python processing.

import sys, os, getopt
import pickle 
input = sys.argv[1:][0]
input  = input + ".pickle"

inputH = open(input) 
data = pickle.load(inputH) 

print "#----END------"
