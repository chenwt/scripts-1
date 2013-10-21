# test, LAML data need to use private key 
# gtdownload -vv -c /Users/jh3283/TCGABAM/pub.key -d 89324e86-5b7a-4f69-92c1-3b67293f8748

# step one is to find the disease 

# search
# cgquery "disease_abbr=LAML" -o LAML.xml
# gtdownload -v -c /Users/jh3283/SCRATCH/school/compGenomic/cghub.key -d ca025fd2-beb4-4327-882f-eca4bb52adcb

cgquery "disease_abbr=LAML&library_strategy=WXS" -o LAML_exome.xml
#366
cgquery "disease_abbr=LAML&library_strategy=RNA-Seq" -o LAML_RNAseq.xml
#180
cgquery "disease_abbr=GBM&library_strategy=RNA-Seq" -o GBM_RNAseq.xml
#339
cgquery "disease_abbr=GBM&library_strategy=WXS" -o GBM_WXS.xml

#1002

cgquery "disease_abbr=LIHC&library_strategy=WXS" -o LIHC_exome.xml
#  219
cgquery "disease_abbr=LIHC&library_strategy=RNA-Seq" -o LIHC_RNAseq.xml
# 57



cgquery "disease_abbr=LAML&library_strategy=WGS" -o LAML_WGS.xml
cgquery "disease_abbr=LIHC&filename=*WGS" -o LIHC_WGS.xml

cgquery "disease_abbr=GBM&library_strategy=WGS" -o GBM_WGS.xml

