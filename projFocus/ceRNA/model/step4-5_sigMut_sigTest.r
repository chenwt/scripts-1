### codes history to calculate fisher exact test of ceRNA driver mutations, ceRNA driver genes for report Jun16_2014
## original data in /ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jun2014/result_step2N3_numbersNtests_06132014.xlsx
## original data in /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene and  /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest

###-----

##### gene level 
#--3utr+2k
fisher.test(x=rbind(c(1391,2329),c(916,1744)))
#--3utr+1k
fisher.test(x=rbind(c(980,2720),c(646,2014)))
#--3utr
fisher.test(x=rbind(c(458,3262),c(299,2361)))
#--2k
fisher.test(x=rbind(c(1095,2625),c(719,1941)))
#--1k
fisher.test(x=rbind(c(624,3096),c(413,2247)))

##### point mutation level
#--1k
fisher.test(x=rbind(c(961,25390),c(613,21090)))
##--2k
fisher.test(x=rbind(c(2014,24337),c(1337,20366)))
##--3utr
fisher.test(x=rbind(c(1072,25279),c(760,20943)))
##--3utr+1k
fisher.test(x=rbind(c(1934,24417),c(1316,20387)))
##--3utr2k
fisher.test(x=rbind(c(2893,23458),c(1983,19720)))

savehistory(file="~/HOME/scripts/projFocus/ceRNA/model/step4-5_sigMut_sigTest.r")

 printFisherTest = function(t) {
   print(c(pval = t$p.value, oddsRatio=t$estimate))
 }
####--------cancer genes ceRNA regulator 
##--gene level
#-- 3utr+2k
t = fisher.test(x=rbind(c(893,1678),c(513,1031)))
printFisherTest(t)
#-- 3utr+1k
t = fisher.test(x=rbind(c(618,1953),c(356,1188)))
printFisherTest(t)
#-- 3utr
t = fisher.test(x=rbind(c(300,2271),c(167,1377)))
printFisherTest(t)
#-- 2k
t = fisher.test(x=rbind(c(691,1880),c(406,1138)))
printFisherTest(t)
#--1k
t = fisher.test(x=rbind(c(378,2193),c(225,1319)))
printFisherTest(t)

##--mutation level
#-- 3utr+2k
t = fisher.test(x=rbind(c(1812,17687),c(1092,13114)))
printFisherTest(t)
#-- 3utr+1k
t = fisher.test(x=rbind(c(1186,18313),c(706,13500)))
printFisherTest(t)
#-- 3utr
t = fisher.test(x=rbind(c(678,18821),c(401,13805)))
printFisherTest(t)
#-- 2k
t = fisher.test(x=rbind(c(1237,18262),c(749,13457)))
printFisherTest(t)
#--1k
t = fisher.test(x=rbind(c(558,18941),c(332,13874)))
printFisherTest(t)
