##------------------------------testing on toy data
setwd("~/Dropbox/compMethod/projects/")
require("plyr")
require("Matrix")
require("MASS")
require("reshape2")
require("useful")

data <- read.table("FEATURE-ICD9.simplified.100patients.ICD9")
colnames(data) <- c("pid","var")
data$value <- 1

dataC <- dcast(data,id ~ var + value,length, value.var="value")

dataMX <- Matrix(as.matrix(dataC),sparse=TRUE)

pdf("ICD9-distrituon.pdf")
image(dataMX, main = "Patients ICD9")
dev.off()

label <- read.table("FEATURE-ICD9.simplified.100patients.pid.LABEL")
colnames(label) <- c("pid","c")