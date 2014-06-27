#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt

## for TF binding site mutation
function getPreTFBS {
  python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step4-1_BSfilter_tf.py
  
  sortxh mut_in_MiRNABindSite_05152014.hg19|uniq > mut_in_MiRNABindSite_05152014.hg19.uniq
  awk 'NR>1{print $1}' mut_in_TFBindSite_05142014.hg19.uniq |sort|uniq > mutGene_in_TFBindStie_05162014
  awk 'NR>1{print $1}' mut_in_MiRNABindSite_05152014.hg19 |sort|uniq > mutGene_in_mirBindSite_05162014
}


###---- all genes
##-- get TFBS mutation in promoter region

function grepMutation {
  fmutIn1k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k
  fmututr3=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p
  ftfbs=$CWD/mut_in_TFBindSite_05142014.hg19.uniq
  fmibs=$CWD/mut_in_MiRNABindSite_05152014.hg19.uniq

  ### $5 field is the mutated gene names
  temp1k=$CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k.temp
  temputr3p=$CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p.temp
  temptf=$CWD/mutGene_in_TFBindStie_05162014.temp
  tempmir=$CWD/mutGene_in_mirBindSite_05162014.temp

  awk '{print $2"-"$6"\t"$0}' $fmutIn1k > $temp1k 
  awk '{print $2"-"$6"\t"$0}' $fmututr3 > $temputr3p 
  awk '{split($2,a,":"); split(a[2],b,"-");print a[1]"-"b[1]}' $ftfbs > $temptf 
  awk '{split($2,a,":"); split(a[2],b,"-");print a[1]"-"b[1]}' $fmibs > $tempmir 
  
  f1ktfbs=$CWD/step5-1_bsSite/allg/ceRNA_driver_greedy_1k.TFBindSite.mut.point_06232014
  futr3pmibs=$CWD/step5-1_bsSite/allg/ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06232014

  ~/bin/grepf2f $temptf $temp1k $f1ktfbs 
  ~/bin/grepf2f $tempmir $temputr3p $futr3pmibs 

   cut -f1 $f1ktfbs |sort|uniq > $f1ktfbs.mutlist 
   cut -f2 $f1ktfbs |sort|uniq > $f1ktfbs.genelist 
   wc -l  $f1ktfbs.mutlist 
   wc -l  $f1ktfbs.genelist 

   cut -f1 $futr3pmibs | sort|uniq > $futr3pmibs.mutlist 
   cut -f6 $futr3pmibs | sort|uniq > $futr3pmibs.genelist 
   wc -l $futr3pmibs.mutlist
   wc -l $futr3pmibs.genelist

   rm $CWD/*temp
}

# grepMutation

###---cancer genes 
##---gene level
function grepMutationByGene {
  fmutIn1k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k
  fmututr3=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p
  
  # awk -F"\t" '{print $5}' $fmutIn1k |sort|uniq > $CWD/cg.1k.mutgene
   ~/bin/grepf2f $CWD/mutGene_in_TFBindStie_05162014  $fmutIn1k $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene_tss.06172014 
   ~/bin/grepf2f $CWD/mutGene_in_mirBindSite_05162014  $fmututr3 $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3.mirBindSite.mutgene_tss.06172014 
  
  ### $5 field is the mutated gene names
  awk '{print $5"\t"$0}' $fmutIn1k > $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k.temp 
  awk '{print $5"\t"$0}' $fmututr3 > $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p.temp 
  
  ~/bin/grepf2f $CWD/mutGene_in_TFBindStie_05162014 $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k.temp  $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014 
   ~/bin/grepf2f $CWD/mutGene_in_mirBindSite_05162014 $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p.temp $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3.mirBindSite.point.mutgene_06172014 
 }
function getMutGene_count {
   cut -f1 $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014 |sort|uniq > $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist
    wc -l $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist 
   cut -f1 $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3.mirBindSite.point.mutgene_06172014 |sort|uniq > $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3.mirBindSite.mutgene_.point.06172014.genelist 
    wc -l $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3.mirBindSite.point.mutgene_06172014.genelist 
}

# getMutGene_count

# grep -wf step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist brca_somaticMutation_all.annotation_Jun172014|awk 'BEGIN{FS=OFS="\t"}{print "chr"$2,$3,$4,$1}' > step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mutgene.point_06172014.genelist.bed


##-- point mutation level

function grepMutation {
  fmutIn1k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k
  fmututr3=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p
  ftfbs=$CWD/mut_in_TFBindSite_05142014.hg19.uniq
  fmibs=$CWD/mut_in_MiRNABindSite_05152014.hg19.uniq

  ### $5 field is the mutated gene names
  temp1k=$CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k.temp
  temputr3p=$CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p.temp
  temptf=$CWD/mutGene_in_TFBindStie_05162014.temp
  tempmir=$CWD/mutGene_in_mirBindSite_05162014.temp

  awk '{print $2"-"$6"\t"$0}' $fmutIn1k > $temp1k 
  awk '{print $2"-"$6"\t"$0}' $fmututr3 > $temputr3p 
  awk '{split($2,a,":"); split(a[2],b,"-");print a[1]"-"b[1]}' $ftfbs > $temptf 
  awk '{split($2,a,":"); split(a[2],b,"-");print a[1]"-"b[1]}' $fmibs > $tempmir 
  
  f1ktfbs=$CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014
  futr3pmibs=$CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014

  ~/bin/grepf2f $temptf $temp1k $f1ktfbs 
  ~/bin/grepf2f $tempmir $temputr3p $futr3pmibs 

   cut -f1 $f1ktfbs |sort|uniq > $f1ktfbs.genelist 
   wc -l  $f1ktfbs.genelist 

   cut -f1 $futr3pmibs | sort|uniq > $futr3pmibs.genelist 
   wc -l $futr3pmibs.genelist

}

# grepMutation

# cd $CWD/step5-1_bsSite
# cat cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014.genelist cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014.genelist |awk -F- '{print "chr"$1"\t"$2"\t"$2+1}' > cg_ceRNA_driver_greedy_1kNutr3.mutation.06172014.bed

# awk '{print $2"-"$3,$1,$4,$5,$6}' brca_somaticMutation_all.annotation_Jun172014 > brca_somaticMutation_all.annotation_Jun172014.temp
# grep -wf $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014.genelist $CWD/brca_somaticMutation_all.annotation_Jun172014.temp | awk -F"-|\t| |;" '
#   BEGIN{OFS="\t"}{print $3,$1,$2,$4,$5,$6,$7,$8,$9}' > $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3pN1k_mitfBS.mutation.anno_Jun172014 
# grep -wf $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014.genelist $CWD/brca_somaticMutation_all.annotation_Jun172014.temp | awk -F"-|\t| |;" '
#   BEGIN{OFS="\t"}{print $3,$1,$2,$4,$5,$6,$7,$8,$9}' >> $CWD/step5-1_bsSite/cg_ceRNA_driver_greedy_utr3pN1k_mitfBS.mutation.anno_Jun172014 

# rm *temp


###-----lasso genes

function grepMutation {
  # fmutIn1k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k
  # fmututr3=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p

  fmutIn1k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k
  fmututr3=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p

  ftfbs=$CWD/mut_in_TFBindSite_05142014.hg19.uniq
  fmibs=$CWD/mut_in_MiRNABindSite_05152014.hg19.uniq

  ### $5 field is the mutated gene names
  temp1k=$CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k.temp
  temputr3p=$CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p.temp
  temptf=$CWD/mutGene_in_TFBindStie_05162014.temp
  tempmir=$CWD/mutGene_in_mirBindSite_05162014.temp

  awk '{print $2"-"$6"\t"$0}' $fmutIn1k > $temp1k 
  awk '{print $2"-"$6"\t"$0}' $fmututr3 > $temputr3p 
  awk '{split($2,a,":"); split(a[2],b,"-");print a[1]"-"b[1]}' $ftfbs > $temptf 
  awk '{split($2,a,":"); split(a[2],b,"-");print a[1]"-"b[1]}' $fmibs > $tempmir 
  
  # f1ktfbs=$CWD/step5-1_bsSite/allg_lasso/ceRNA_driver_lasso_1k.TFBindSite.mut.point_06232014
  # futr3pmibs=$CWD/step5-1_bsSite/allg_lasso/ceRNA_driver_lasso_utr3p.mirBindSite.mut.point_06232014

  f1ktfbs=$CWD/step5-1_bsSite/cg_lasso/ceRNA_driver_lasso_1k.TFBindSite.mut.point_06232014
  futr3pmibs=$CWD/step5-1_bsSite/cg_lasso/ceRNA_driver_lasso_utr3p.mirBindSite.mut.point_06232014


  ~/bin/grepf2f $temptf $temp1k $f1ktfbs 
  ~/bin/grepf2f $tempmir $temputr3p $futr3pmibs 

   cut -f1 $f1ktfbs |sort|uniq > $f1ktfbs.mutlist 
   cut -f2 $f1ktfbs |sort|uniq > $f1ktfbs.genelist 
   wc -l  $f1ktfbs.mutlist 
   wc -l  $f1ktfbs.genelist 

   cut -f1 $futr3pmibs | sort|uniq > $futr3pmibs.mutlist 
   cut -f6 $futr3pmibs | sort|uniq > $futr3pmibs.genelist 
   wc -l $futr3pmibs.mutlist
   wc -l $futr3pmibs.genelist

   rm $CWD/*temp
}

grepMutation


