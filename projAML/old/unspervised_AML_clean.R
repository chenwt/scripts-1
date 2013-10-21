setwd("Y:/AML")

exp.data <- read.table("target_AML_uniq_probes_177.exp", header=T)
setwd("C:/Documents and Settings/jh3283/Desktop/AML/")
save.image("target_uniq_177_exp.Rdata")

require("ggplot2")

data <- na.omit(exp.data)
rownames(data) <- data[,2]
names(data)
data<- data[,-c(1:2)]

##### cluster on samples
exp.t <- t(data)
dim(exp.t)
row.names(exp.t)
x <- exp.t

#####
# Ward Hierarchical Clustering
x <- na.omit(x)
d <- dist(x, method = "euclidean") 
fit2 <- hclust(d, method="complete")
groups <- cutree(fit2, k=3) 

pdf("cluster.pdf")
plot(fit2,cex.axis=0.5,labels=F,xlab = "groups") 
rect.hclust(fit2, k=3, border="blue")
dev.off()

ls.groups <- as.list(groups)
group1 <- names(ls.groups[which(ls.groups=="1")])
group2 <- names(ls.groups[which(ls.groups=="2")])
group3 <- names(ls.groups[which(ls.groups=="3")])

exp <- exp.data
exp.g1 <- cbind(exp[,c(2,1)],exp[,group1])
exp.g2 <- cbind(exp[,c(2,1)],exp[,group2])
exp.g3 <- cbind(exp[,c(2,1)],exp[,group3])


write.table(exp.g1,"AML_group1_68.exp",sep="\t",row.names = F )
write.table(exp.g2,"AML_group2_23.exp",sep="\t")
write.table(exp.g3,"AML_group3_86.exp",sep="\t")

groupsDF <- as.data.frame(unlist(ls.groups))

write.table(groupsDF,"sample_group.txt",col.names = F)
