find exome/ -type f -ls |grep bam |awk 'BEGIN{OFS="\t";} {print $11,$7}' > CHECK/files_GBM_exome_downloaded.txt
awk 'BEGIN{OFS="\t";}{print "exome/"$1"/"$2,$3}' CHECK/GBM_exome.txt > CHECK/files_GBM_exome_want.txt
grep --file=files_GBM_exome_downloaded.txt files_GBM_exome_want.txt  > files_GBM_exome_sucessDownloaded.txt

grep --file=temp.txt GBM_exome.txt |awk 'BEGIN{OFS="\t";}{if($4==01&&$5=="HG19_Broad_variant") print "exome/"$1"/"$2}'|sort -nk 3 > files_GBM_exome_sucessTumor_HG19.txt

grep --file=temp.txt GBM_exome.txt |awk 'BEGIN{OFS="\t";}{if($4==10&&$5=="HG19_Broad_variant") print "exome/"$1"/"$2}'|sort -nk 3  > files_GBM_exome_sucessNormal_HG19.txt
cat files_GBM_exome_sucessTumor_HG19.txt | cut -d/ -f3 | cut -d'-' -f3 | sort > pIDTumor.temp
cat files_GBM_exome_sucessNormal_HG19.txt | cut -d/ -f3 | cut -d'-' -f3 | sort > pIDNormal.temp
grep --file=pIDNormal.temp pIDTumor.temp > pairID.temp
grep --file=pairID.temp files_GBM_exome_sucessTumor_HG19.txt > files_GBM_exome_sucess_tumor_pair.txt
grep --file=pairID.temp files_GBM_exome_sucessNormal_HG19.txt > files_GBM_exome_sucess_normal_pair.txt
