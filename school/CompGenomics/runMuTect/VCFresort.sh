cp $1 $1.temp
cat $1 | grep "^#" > header.temp 
sed '/^#/d' $1.temp > $1.temp.temp
# for $CHR in {1..9 X Y 10..22}
# do 
# 	cat $1.temp.temp | awk -v awk_chr=$CHR '{if($1==awk_chr) print $0}'  >> $1.resorted
# done

cat $1.temp.temp | awk '{if($1>=1 && $1<10) print $0}'  > $1.temp

sort -n -k 1,1 -n -k 2,2 < $1.temp > $1.temp.1
sed '/^[1-9]\t/d' $1.temp.temp  > $1.temp

cat $1.temp.temp | awk '{if($1=="X") print $0}'  > $1.temp
sort -k 1,1 -k 2,2 < $1.temp > $1.temp.2

cat $1.temp.temp| awk '{if($1=="Y") print $0}'  > $1.temp
sort -n -k 1,1 -n -k 2,2 < $1.temp > $1.temp.3

cat $1.temp.temp | awk '{if($1>9 && $1<23) print $0}'  >  $1.temp
sort -n -k 1,1 -n -k 2,2 < $1.temp > $1.temp.4

# sed '/^[12XY]/d' $1.temp.temp  > $1.temp

cat header.temp $1.temp.1 $1.temp.2 $1.temp.3 $1.temp.4 > $1.temp
uniq -u $1.temp > $1.resorted
# still need to remove those lines which have the same ref and alt
rm *temp*