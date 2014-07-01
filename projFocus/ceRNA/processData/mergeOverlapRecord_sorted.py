#!/usr/bin/python
#J.HE
'''
input: sorted bed file
output: bed file with overlapped record merges, merged recoded countted
purpose: data process
status: underwork

'''

import os, sys, re

argv = sys.argv[1:] 
if len(argv) == 1:
    output = input + ".merged"
else: 
    input, output = argv[:2]

fhout = open( output, 'w') 

reg_prev = ['0', 0, 0 ] 
with(open(input)) as f:
    line = f.readline()
    while line:
        chr, ps, pe, info = line.strip().split("\t",3)
        if chr = reg_prev[0]:

        line = f.readline()

