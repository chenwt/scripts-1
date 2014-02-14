#!/bin/bash
#By: J.He
#$ -cwd 

cnt=0
while [ $cnt -lt 3 ]; do
  let cnt=$cnt+1
  echo "$cnt test" > file$cnt.txt
  echo $cnt
  sleep 2s 
done

echo $cnt
sleep 2s
cntNF=`ls -1 file*txt 2>/dev/null | wc -l`
cntjob=0
echo "new file.." $cntNF
while [ $cntNF -lt $cnt ] 
      [ $cntjob -lt $cnt ]
do
    if [ $cntNF -gt 1 ] ; then
      for line in `ls file*txt` 
      do
          cat $line > secondJobFile_$line
	  ((cntjob += 1))
          echo "new file generated $cntjob " >> secondJobFile_$line
          rm $line 
          echo "new file generated $cntjob"
      done
    elif [ $cntNF -eq 1 ]; then 
      line=`ls flle*txt`
      cat $line > secondJobFile_$line
      ((cntjob += 1))
      echo "new file generated $cntjob " >> secondJobFile_$line
      rm $line
      echo "new file generated $cntjob"
    else
      echo "waiting..."
    fi
    sleep 3s
    cntNF=`ls -1 file*txt 2>/dev/null | wc -l`
done

echo $cntjob"new file generated"
echo "#---DONE"
