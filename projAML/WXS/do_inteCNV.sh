#input the subfolder name of resTN or resRN
#get the information of specific PID
sd=~/script/projAML/WXS/
wd=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/CNV/
od=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/CNV/report/
fname="ALL.DATA.cnv"

cd ${wd}
head -1 $1/${fname} > out.temp
while read pid 
do
	grep ${pid} $1/${fname} >> out.temp
done < ${pidf}