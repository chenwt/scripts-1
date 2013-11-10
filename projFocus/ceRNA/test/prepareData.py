#!/usr/bin/python
#input: <string: gz file name from TCGA>
#output: NO
#Desp: 	1. unzip the gz file and delete uncessary ones; 
# 		2. rename all files using TCAG barcode(from the sdrf file)
# 		3. generate inputfile for getMat.py()
# 		4. run getMat.py to get the big matrix
# 		5. remove all uncessary files
