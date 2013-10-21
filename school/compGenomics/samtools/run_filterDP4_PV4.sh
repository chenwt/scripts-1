cat ../tumor.Y.flt.vcf | grep -v ^# |cut -f8 | sed -e 's/\;/\t/g' | cut -f4 | sed -e 's/DP4\=//g' |sed -e 's/\,/\t/g' | awk '{if(>2&&>2&&+>4) print -sh}' | head
