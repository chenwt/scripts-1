#SCRIPT:  school/compGenomics/gtBatch.sh
#PURPOSE: Process a file line by line with redirected while-read loop.
#need to change the key address if want to use it in local

FILENAME=$1
count=0
KEY="/ifs/scratch/c2b2/ac_lab/jh3283/school/compGenomic/key/cghub.key"
while read LINE
do	
	let count++
      	echo "$count $LINE" >> logs.downloaded
      	cmd="gtdownload -v -c $KEY -d $LINE"
	echo $cmd >> logs.download 
	$cmd >> logs.download 
	# echo "$LINE" >> logs.downloaded
done < $FILENAME
echo -e "\nTotal $count lines read"
