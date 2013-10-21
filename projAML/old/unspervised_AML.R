setwd("Y:/AML")

exp.data <- read.table("target_AML_uniq_probes_177.exp", header=T)

setwd("C:/Documents and Settings/jh3283/Desktop/AML/")
save.image("target_uniq_177_exp.Rdata")

source("http://bioconductor.org/biocLite.R")
biocLite("affy")
library("affy")
install.packages("RColorBrewer")
install.packages("pheatmap")
require("RColorBrewer")
require("pheatmap")
require("ggplot2")

plot.density(x)

data <- na.omit(exp.data)
rownames(data) <- data[,1]
names(data)
data<- data[,-c(1:2)]
data <- apply(data,2,FUN=as.numeric)
summary(data[,1:10])

heatmap(data)
##### cluster on samples

exp.t <- t(exp.data)
dim(exp.t)
row.names(exp.t)



x <- exp.t

row.names(x) <- gsub(".CEL","",row.names(x))
row.names(x) <- gsub("X2010.06.25","X0625",row.names(x))
row.names(x) <- gsub("X20100623","X0623",row.names(x))
row.names(x) <- gsub("1037","",row.names(x))
row.names(x) <- gsub("_",".",row.names(x))

colnames(x)
colnames(x) <- x[1,]
x <- x[2:length(x[,1]),]
x <- apply(x,c(1,2),FUN=as.numeric)
x <- na.omit(x)
length(x[,1])
dim(x)
x[1:5,1:5]


#####  k value testing
k_max <- 30
eval_k <- function(x,k_max)
{
	wss <- data.frame(k=c(1:k_max),wss=(nrow(x)-1)*sum(apply(x,2,var)),flag=1)
    for (i in 2:k_max) wss$wss[i] <- sum(kmeans(x, centers=i)$withinss)
	for (i in 2:k_max) wss$flag[i] <- (wss$wss[i]-wss$wss[i-1])/wss$wss[i-1]
	return(wss)
}
# K-Means Cluster Analysis
myplot <- function(wss,k_max){
	p1 <- ggplot(wss,aes(x= 1:k_max,y=wss)) + geom_line(colour = "lightblue",size = 1) + geom_point(pch = 20, colour = "blue",size =3) + xlab("k") + opts(title="testing K")
	ik <- NA
    iwss <-NA
    ik<- wss$k[which(wss$flag>0)]
	iwss<- wss$wss[which(wss$flag>0)] 
    
		if(length(ik)>0)
		{
			p1 <- p1 + geom_point(aes(x=ik, y=iwss,color="red", pch = 20,size =3))
		}
	return(p1)
}

p1

fit <- kmeans(x, 11) 

install.packages("cluster")
library(cluster) 
x[1:5,1:5]
dim(x)
clusplot(pam(x,11),col.p="blue")
# clusplot(x,fit$cluster, color=TRUE, shade=TRUE, labels=11, lines=0)
install.packages("fpc")
library(fpc)
plotcluster(x, fit$cluster)


aggregate(x,by=list(fit$cluster),FUN=mean)
x <- data.frame(x, fit$cluster)
ord<- order(x[,dim(x)[2]])
x <- x[ord,]
pheatmap(x[,-dim(x)[2]],cluster_rows=F,cluster_cols=F,col=brewer.pal(11,"Set3"),border_color = NA)
dim(x)


#####
# Ward Hierarchical Clustering
x <- na.omit(x)
d <- dist(x, method = "euclidean") 
fit <- hclust(d, method="ward") 
fit2 <- hclust(d, method="complete")
fit3 <- hclust(d,method="average")
plot(fit2,cex.axis=0.5,labels=F,xlab = "groups") 
groups <- cutree(fit2, k=3) 
rect.hclust(fit2, k=3, border="blue")
fit$order
attributes(fit)
fit2$order

ls.groups <- as.list(groups)
group1 <- names(ls.groups[which(ls.groups=="1")])
group2 <- names(ls.groups[which(ls.groups=="2")])
group3 <- names(ls.groups[which(ls.groups=="3")])

#######################

# DEG

##########
exp <- t(x)
sample <- colnames(exp)

exp.g1 <- exp[,group1]
exp.g2 <- exp[,group2]
exp.g3 <- exp[,group3]

var1 <- var(exp.g1)
var2 <- var(exp.g2)
var3 <- var(exp.g3)

res <-   data.frame(t12=rep(0,length(exp.g1[,1])),row.names = rownames(exp.g1))

res <- list()
Calt <- function(exp.g1,exp.g2)
{
 for (i in 1: length(exp.g1[,1]))
    {
     	temp <- t.test(exp.g1[i,],exp.g2[i,])
     	res[[i]] <- as.numeric(temp$statistic)
     }
    return (res)
}

result.t12 <- Calt(exp.g1,exp.g2)