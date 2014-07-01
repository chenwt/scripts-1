CWD = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/bs_with_addedBS"

fgred = paste(CWD, "/", "ceRNA_driver_greedy_mut.mat.uniq.bed", sep= "")
flass = paste(CWD, "/", "ceRNA_driver_lasso_mut.mat.uniq.bed", sep ="")

gredD = read.table(fgred, sep="\t")
lassD = read.table(flass, sep="\t")
gredD$V1 <- vapply(gsub("X", "23", gredD$V1),as.integer, 1)
lassD$V1 <- vapply(gsub("X", "23", lassD$V1), as.integer, 1)

for ( i in c(2,3)) {
  gredD[,i] <- vapply(gredD[,i],as.numeric, 1)
  lassD[,i] <- vapply(lassD[,i],as.numeric, 1)
}

plot(y = lassD$V1, x = lassD$V2, col="black", type = "p", pch = 16)
points(y = gredD$V1, x = gredD$V2, col="red", pch = 16, cex = 0.6)

barplot(horiz=T)
