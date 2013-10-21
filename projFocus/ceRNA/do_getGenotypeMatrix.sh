#!/bin/bash
#By: J.He
#TODO: 
#input: <full path of filenames of TCGA level 2 SNP birdseed data>
#output:<file of snp * sample>

# datad=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/
resd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result_snp/
fileout="testout.birdseed"
filein="test.birdseed"
# echo $datad
echo $resd
echo $fileout
# echo $input
echo $1
# mkfifo pipein pipeout
cnt=1
temp="tempFile"
one=1
while read filein
do 
	if [[ ${cnt} -eq 1 ]]; then
	# head ${datad}${filein}
		echo "first file"
		awk 'BEGIN{
				OFS=FS="\t"
			}
			{if(NR < 2) 
			 	print $1,$2
			else if (NR > 2)
			  	if($3 > 0.01)
				 	print $1,"-1"
				else
					print $1,$2
			}' ${filein} > ${resd}${fileout}
		head ${resd}${fileout}
		cnt=$(echo ${cnt}+1|bc)
	else 
		# echo $cnt"File"
		awk 'BEGIN{
				OFS=FS="\t"
			}
			{if(NR < 2) 
			 	print $2
			else if (NR > 2)
			  	if($3 > 0.01)
				 	print "-1"
				else
					print $2
			}' ${filein} > ${resd}${temp}
		head ${resd}${temp}
		paste ${resd}${fileout} ${resd}${temp} > ${restd}${temp}2
		mv ${restd}${temp}2 ${resd}${fileout}
		head ${resd}${fileout}
		# cat pipeout > ${resd}${fileout}
		cnt=$(echo ${cnt}+1|bc)
		# echo $cnt
	fi
	echo $cnt
done < $1

# rm pipein pipeout
## }' ${datad}${filein} > ${restd}${fileout}
	# head ${restd}${fileout}
	