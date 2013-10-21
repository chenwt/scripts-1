require(useful)

setwd("Z:././AML/level3/")

# loading expression data(normalized)
data <- read.delim("TARGET_RMA_norm_GE_Level3.txt",header=TRUE)
colnames(data) <- tolower(colnames(data))
data <- data[,-c(1,3)]
# loading centroids data from valk2 16 cluster, format gene name + value
centroids <- read.csv("../valk2/centroids.shrunk.csv",header=T)
colnames(centroids) <- tolower(colnames(centroids))
centroids <- centroids[,-1]
length(unique(data$level))==nrow(data)
#centroid.overall <- read.csv("../valk2/centroids.overall.csv",header=T)
#colnames(centroid.overall) <- tolower(colnames(centroid.overall))
#centroid.overall <- centroid.overall[,-2]
geneSet <- unique(centroids$annotation)
#apply(genes,)

# generate the intersection
 gene.intersect <- intersect(geneSet,data$level)
		#test
		tmp <- gene.intersect[1:10]
      
		tmp2 <- centroids[which(centroids$annotation == tmp[1]),]
		tmp2[names(which(apply(tmp2[,3:ncol(tmp2)],1,mean) == max(apply(tmp2[,3:ncol(tmp2)],1,mean)))),]

 uniqgene_centroid <- function(centroids,geneset){
 		res <- data.frame()
		for ( i in 1:length(geneset)) {
		    ## gene names in the first colum
	        tmp <- centroids[which(centroids[,1] == geneset[i]),]
            if(is.null(dim(centroids[,-1]))) {
                res <-rbind(res,tmp[names(which((tmp[,-1]) == max(tmp[,-1]) ) ),])
            }else{
                res <-rbind(res,tmp[names(which(apply(tmp[,-1],1,mean) == max(apply(tmp[,-1],1,mean)))),])
	        }
		}
		return (res) 	
 }
centroids.new <- uniqgene_centroid(centroids,geneSet)

# clustering level3 data using 16 centroids

 test


# save data
save.image("level3_pamr.RData")