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
rownames(data) <- data[,2]
data<- data[,-c(1:2)]
data <- apply(data,2,FUN=as.numeric)
summary(data[,1:10])

heatmap(data)
##### cluster on samples
exp.t <- t(exp.data)
dim(exp.t)
row.names(exp.t)



x <- exp.t
x <- na.omit(x)
names(x)
x[1:3,1:5]
k_max <- 30

#####  k value testing
eval_k <- function(x,k_max)
{
	wss <- data.frame(k=c(1:k_max),wss=(nrow(x)-1)*sum(apply(x,2,var)),flag=1)
    for (i in 2:k_max) wss$wss[i] <- sum(kmeans(x, centers=i)$withinss)
	for (i in 2:k_max) wss$flag[i] <- (wss$wss[i]-wss$wss[i-1])/wss$wss[i-1]
	return(wss)
}
# K-Means Cluster Analysis
fit <- kmeans(x, 11) 

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



aggregate(x,by=list(fit$cluster),FUN=mean)
x <- data.frame(x, fit$cluster)
ord<- order(x[,dim(x)[2]])
x <- x[ord,]
pheatmap(x[,-dim(x)[2]],cluster_rows=F,cluster_cols=F,col=brewer.pal(11,"Set3"),border_color = NA)
dim(x)


#####
# Ward Hierarchical Clustering
d <- dist(x, method = "euclidean") 
fit <- hclust(d, method="ward") 
plot(fit) 
groups <- cutree(fit, k=5) 
rect.hclust(fit, k=5, border="red")

# Ward Hierarchical Clustering with Bootstrapped p values
tmp.x<- sample(x)

library(pvclust)
fit <- pvclust(x, method.hclust="ward", method.dist="euclidean")
plot(fit) 
pvrect(fit, alpha=.95)


library(mclust)
fit <- Mclust(x)
plot(fit, x) # plot results 
print(fit) # display the best model

# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 11)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster) 
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
  	 labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mydata, fit$cluster)