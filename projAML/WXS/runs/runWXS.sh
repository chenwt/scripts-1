#!/bin/bash
#By: J.He
#TODO: 


grep -i -f PID_16.txt ~ys2559/data/AML/analysis/predictedRisk.txt|cut -d' ' -f2|sort|uniq -c

