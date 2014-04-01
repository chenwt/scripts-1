#!/bin/bash
#! -cwd
#By: J.He
#Desp.: given a genelist, return the associated network in given network database
#


#----given a genename list, get preppi network 
# Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/network/geneName2Uniprot.r publicGene.list
# Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/network/geneName2Uniprot.r yLanGene.list

# source /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/network/plotPPI.r 


#------
# CWD=/ifs/home/c2b2/ac_lab/jh3283/DATA/projMisc/yLan/Mar4
SRCDIR=/ifs/home/c2b2/ac_lab/jh3283/scripts/
PREPPIDB=/ifs/scratch/c2b2/ac_lab/jh3283/database/preppi/preppi_int_3col_genename.txt_600_90
END="[#END---]"
cd $CWD

preppiNet(){
    genelist=$1
    grep -w -f $genelist $PREPPIDB > $genelist.preppi
    ~/tools/R/R_current/bin/Rscript $SRCDIR/projAML/WXS/network/plotPPI.r $genelist 
    echo $END
}

# genelist=$CWD/wesGene.list
# preppiNet $genelist 
# genelist=$CWD/combinded_WesGene_PublicGene.list
# preppiNet $genelist 

##-----Mar 29 2014
CWD=/ifs/home/c2b2/ac_lab/jh3283/DATA/projMisc/yLan/Mar29
preppiNetfrom2list(){
    glist1=$1
    glist2=$2
    # grep -w -f $glist1 $PREPPIDB > $glist1.preppi
    # grep -w -f $glist2 $PREPPIDB > $glist2.preppi
    ~/tools/R/R_current/bin/Rscript $SRCDIR/projMisc/yLan/plotPPI2list.R $glist1 $glist2 
    echo $END
}

l1=$CWD/denovalGenename.list
l2=$CWD/mouseGenename.list
preppiNetfrom2list $l1 $l2 

l2=$CWD/mouseHPgenename.list
preppiNetfrom2list $l1 $l2 
l2=$CWD/mouseTPgenename.list
preppiNetfrom2list $l1 $l2 
