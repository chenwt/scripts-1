###


pdf("1.pdf",width=5)
plot.table(regulators,text.cex=5)
dev.off()

##-----finding class intervals for continious numerical variables
install.packages("classInt")

