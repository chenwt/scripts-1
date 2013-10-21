
awk 'BEGIN{OFS="\t"} {print $1"_"$2,$3,$4,FILENAME}' PASFEW-*.freq.temp |sort -k1 
awk 'BEGIN{OFS="\t"} {print $1"_"$2,$3,$4,FILENAME}' PASFEW-TuA.freq.temp 
awk 'BEGIN{OFS="\t"} {print $1"_"$2,$3,$4,FILENAME}' PASFEW-ReA.freq.temp 