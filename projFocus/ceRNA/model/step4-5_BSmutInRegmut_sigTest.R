### codes history to calculate fisher exact test of ceRNA driver mutations, ceRNA driver genes for report Jun16_2014
## original data in /ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jun2014/result_step5-1_functionalFiltering_Jun2014.xlsx

###-----all genes
##### point mutation level
 printFisherTest = function(t) {
   oddsRatio = ifelse(t$estimate < 1, 1/t$estimate, t$estimate)
   print(c(pval = t$p.value, oddsRatio=oddsRatio))
 }
#--1k TFBS
t = fisher.test(x=rbind(c(23980,19955),c(1001,702)))
printFisherTest(t)
##--3utr mirBS
 t = fisher.test(x=rbind(c(23980,19955),c(1504,1133)))
 printFisherTest(t)
##--3utr+1k
 t = fisher.test(x=rbind(c(23980,19955),c(2370,1747)))
 printFisherTest(t)


##--
#--1k TFBS
t = fisher.test(x=rbind(c(2118,1562),c(111,74)))
printFisherTest(t)
##--3utr mirBS
 t = fisher.test(x=rbind(c(2118,1562),c(142,112)))
 printFisherTest(t)
##--3utr+1k
 t = fisher.test(x=rbind(c(2118,1562),c(252,185)))
 printFisherTest(t)

 
 #### binding site
 t = fisher.test(x=rbind(c(25275,20852),c(244,179)))
 printFisherTest(t)
 ##--3utr mirBS
 t = fisher.test(x=rbind(c(25275,20852),c(872,707)))
 printFisherTest(t)
 ##--3utr+1k
 t = fisher.test(x=rbind(c(25275,20852),c(1075,850)))
 printFisherTest(t)

 
 #### binding site
 t = fisher.test(x=rbind(c(17965,13196),c(651,402)))
 printFisherTest(t)
 ##--3utr mirBS
 t = fisher.test(x=rbind(c(17965,13196),c(981,667)))
 printFisherTest(t)
 ##--3utr+1k
 t = fisher.test(x=rbind(c(17965,13196),c(1533,1009)))
 printFisherTest(t)
 
 
 t = fisher.test(x=rbind(c(17965,13196),c(651,402)))
 printFisherTest(t)
 ##--3utr mirBS
 t = fisher.test(x=rbind(c(17965,13196),c(981,667)))
 printFisherTest(t)
 ##--3utr+1k
 t = fisher.test(x=rbind(c(17965,13196),c(1533,1009)))
 printFisherTest(t)
 
 
 t = fisher.test(x=rbind(c(18734,13658),c(157,102)))
 printFisherTest(t)
 ##--3utr mirBS
 t = fisher.test(x=rbind(c(17965,13196),c(630,458)))
 printFisherTest(t)
 ##--3utr+1k
 t = fisher.test(x=rbind(c(17965,13196),c(764,547)))
 printFisherTest(t)
