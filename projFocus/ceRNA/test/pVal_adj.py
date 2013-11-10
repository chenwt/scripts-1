#!/usr/bin/python
#input: <file: output from Test_snp_exp.py>
#output:<file: with adjusted p-value, and sorted by the adjusted p-value>
#TODO: try to use pickle.

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import numpy as np
stats = importr('stats')
fname = 'test.txt'
data = np.genfromtxt(fname, delimiter=',', names=True)
pvalue_list = data[:,4] #5th column is always p-value
p_adjust = stats.p_adjust(FloatVector(pvalue_list), method = 'BH')

