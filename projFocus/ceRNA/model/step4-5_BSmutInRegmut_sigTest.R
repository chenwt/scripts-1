### codes history to calculate fisher exact test of ceRNA driver mutations, ceRNA driver genes for report Jun16_2014
## original data in /ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jun2014/result_step5-1_functionalFiltering_Jun2014.xlsx

###-----all genes
##### point mutation level
 printFisherTest = function(t) {
   oddsRatio = ifelse(t$estimate < 1, 1/t$estimate, t$estimate)
   print(c(pval = t$p.value, oddsRatio=oddsRatio))
 }
#--1k TFBS
t = fisher.test(x=rbind(c(92,869),c(57,556)))
printFisherTest(t)
##--3utr mirBS
t = fisher.test(x=rbind(c(69,1003),c(64,696)))
printFisherTest(t)
##--3utr+1k
t = fisher.test(x=rbind(c(161,1773),c(121,1195)))
printFisherTest(t)


##--cancer gene 
#--1k TFBS
t = fisher.test(x=rbind(c(55,503),c(31,301)))
printFisherTest(t)
##--3utr mirBS
t = fisher.test(x=rbind(c(46,632),c(34,367)))
printFisherTest(t)
##--3utr+1k
t = fisher.test(x=rbind(c(101,1085),c(65,641)))
printFisherTest(t)

