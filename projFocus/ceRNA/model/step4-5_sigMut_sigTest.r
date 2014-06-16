### codes history to calculate fisher exact test of ceRNA driver mutations, ceRNA driver genes for report Jun16_2014


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
