#!/usr/bin/python
#J.HE
'''
Desp.: pipeline to run greedy, 
'''


import sys,getopt
argv = sys.argv[1:]

usage = 'python ' + sys.argv[0] + ' \n \
        -L <level=, runnning level: \
            0 for data checking; \
            1 for gene-samples list;\
            2 for group lasso;\
            3 for greedy>\n  \
        -I <inputG=, gene to study, in a list file, must be in network file> \n\
        -D <dataDir=, data folder including:\
            1. network data with 1,2 as regulaor-target;\
            2. expression data, first row is sample name, first column is gene symbolei\
            3. somatic mutation data  \
            4. cnv matrix (optional)  \
            5. methylation (optional)  \n\
        -O <outDir=, output dir>\n \
        -T <tempDir=, temp file dir,(optional), if no, use outputd > \n\
        '

example = 'python ' + sys.argv[0] + ' \n \
        -N <netType=, network type: \
            cerna, no direction, no negative correlation; \
            tf, have direction, both negative or positive \
            tfmd, tf-modulator network \n  \
        -I <inputG=, gene to study, in a list file, must be in network file> \n\
        -D <dataDir=, data folder including:\
            1. network data with 1,2 as regulaor-target;\
            2. expression data, first row is sample name, first column is gene symbolei\
            3. somatic mutation data  \
            4. cnv matrix (optional)  \
            5. methylation (optional)  \n\
        -O <outDir=, output Dir>\n \
        -T <tempDir=, temp file dir,(optional), if no, use outputd > \n\
        '

##--get arguments
try:
    opts,args = getopt.getopt(argv,"hI:N:D:O:T:")
except getopt.GetoptError:
    print usage + "\n" + example 
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print usage + "\n" + example 
        sys.exit()
    elif opt in ("-N"):
        networkType = arg
    elif opt in ("-D"):
        dataDir = arg
    elif opt in ("-I"):
        tgeneFile = arg
    elif opt in ("-O"):
        outputDir = arg
    elif opt in ("-T"):
        tempDir = arg


##--QC : check file, folder legitamation

##--step 1: get gslist 

##--step 2: run regreesion 

##--step 3: run greedy  

##--step 4: generate reports
