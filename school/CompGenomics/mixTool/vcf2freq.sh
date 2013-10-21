
PID=$1
cd $PID

if [ -f $PID.somaticfiltered.freq ] ;then  rm $PID.somaticfiltered.freq ; fi
for chrom in `seq 1 22` X Y  
do
	cat $chrom.var.vcf.somaticfiltered.vcf | grep -v "^#"| perl -lne '@F=split("\t");$F[7]=~m/DP4=(.*?);/;$PF=$1; $PF =~ s/,/\t/g;print "$F[0]\t$F[1]\t$PF"' | awk 'BEGIN{OFS="\t"}{print $1,$2,$3+$4,$5+$6}' >> ../$PID.somaticfiltered.freq
done
echo $PID"freq file created!" >>logs$PID


