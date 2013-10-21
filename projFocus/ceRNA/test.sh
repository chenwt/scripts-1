#!/bin/bash
#$-cwd
# J.HE
wd=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/test/
filein=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/test1.txt
fileout=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/test_joint.txt

mkfifo pipein pipeout
for cnt in $(seq 1 3)
do 
	if [ ${cnt}=='1' ]; then
		
	${wd}? > pipein
