# TODO: Add comment
# 
# Author: fmg2117
###############################################################################


stouffer<-function(x){
	Z  <- sum(x)/sqrt(length(x))
	return(Z)
}

z2p<-function(z){
	pnorm(abs(z), lower.tail=F)*2
}



loopClustering<-function(
		vipermat=NULL,
		mut_pvalue=NULL,
		mut_assoc=NULL,
		entrez2gene=NULL,
		tumorexpmat=NULL,
		normexpmat=NULL,
		regul=NULL,
		maxstep=20,
		maxmrs=20,
		maxmutations=20,
		activityThreshold=1.96,
		ledge=FALSE # should the samples included be the leading edge ones?
) {
	
	
	
	### Initialize the loop to create cluster seeds
	
	
	# These objects will be reduced as we build new clusters
	mut_pvalue_here<-mut_pvalue
	mut_assoc_here<-mut_assoc
	vipermat_here<-vipermat
	
	# Summarization object containing a list of the subclusters
	clusterlist<-list()
	
	
	
	step<-0
	skippedStep<-FALSE # For failsafe skip of a step
	while(TRUE){
		if(skippedStep){
			skippedStep<-FALSE			
		} else{
			step<-step+1
		}
		if(step>maxstep){
			return(clusterlist)
		}
		
		# Starting with the best GenomicFeature-activity pairs
		best_muts<-which(mut_pvalue_here == min(mut_pvalue_here), arr.ind = TRUE)
		
		# Best cluster
		seed<-rownames(best_muts)[1]
		cat("Now processing ",entrez2gene[seed],"\n",sep="")
		seedp<-mut_pvalue_here[best_muts[1,1],best_muts[1,2]]
		if(seedp>0.05){
			cat(paste("All samples were assigned to a mutation/activity cluster at step ",step,"\n",sep=""))
			return(clusterlist)
		}
		
		
		# Select the mutations significantly associated to this gene (or the top maxmutations sorted by pvalue)
		these_muts<-mut_pvalue_here[best_muts[seed,1],]
		these_muts<-these_muts[these_muts<=0.01]
		if(length(these_muts)>=maxmutations) {
			top_muts<-sort(these_muts)[1:maxmutations]
		} else {
			top_muts<-these_muts
		}
		
		# Generate column association (mutation) color object
		colAnnotation<-matrix("white",nrow=length(top_muts),ncol=ncol(vipermat_here))
		rownames(colAnnotation)<-names(top_muts)
		colnames(colAnnotation)<-colnames(vipermat_here)
		for (mutgene in names(top_muts)) {
			origene<-mutgene
			if(length(grep("^SNV_",mutgene,perl=TRUE))==1){
				mutgene<-sub("^SNV_","",mutgene,perl=TRUE)
				col<-rep("white",length(rawsnp[mutgene,]))
				col[rawsnp[mutgene,]!=0]<-"black"
				names(col)<-names(rawsnp[mutgene,])
			} else if(length(grep("^AMP_",mutgene,perl=TRUE))==1){
				mutgene<-sub("^AMP_","",mutgene,perl=TRUE)
				col<-rep("white",length(rawfcnv[[mutgene]]))
				col[rawfcnv[[mutgene]]>0]<-"red"
				col[rawfcnv[[mutgene]]<0]<-"blue" # TODO check why this is needed? # TODO combine them?
				names(col)<-names(rawfcnv[[mutgene]])
			} else if(length(grep("^DEL_",mutgene,perl=TRUE))==1){
				mutgene<-sub("^DEL_","",mutgene,perl=TRUE)
				col<-rep("white",length(rawfcnv[[mutgene]]))
				col[rawfcnv[[mutgene]]>0]<-"red"
				col[rawfcnv[[mutgene]]<0]<-"blue"
				names(col)<-names(rawfcnv[[mutgene]])
			} else if(length(grep("^MET_",mutgene,perl=TRUE))==1){
				mutgene<-sub("^MET_","",mutgene,perl=TRUE)
				col<-rep("white",length(rawfmeth[[mutgene]]))
				col[rawfmeth[[mutgene]]!=0]<-"green"
				names(col)<-names(rawfmeth[[mutgene]])
			}
			commonsamples<-intersect(names(col),colnames(colAnnotation))
			colAnnotation[origene,commonsamples]<-col[commonsamples]
		}
		rownames(colAnnotation)<-paste(sep="",substr(rownames(colAnnotation),1,4),entrez2gene[substr(rownames(colAnnotation),5,nchar(rownames(colAnnotation)))])
		
		
		# Keep only samples harboring a mutation in our top list
		mutationSamples<-c()
		for(j in 1:ncol(colAnnotation)) {
			if(all(colAnnotation[,j]=="white")){
				
			} else {
				mutationSamples<-c(mutationSamples,colnames(colAnnotation)[j])
			}
		}
		
		# For this gene, get only samples within the leading edge of enrichment between the seed and its most significant mutation
		if(ledge){
			# Get the mutated samples for the top mutation
			top_mut<-names(top_muts)[1]
			if(length(grep("SNV_",top_mut))==1){
				top_mut_symbol<-sub("SNV_","",top_mut)
			}else if(length(grep("AMP_",top_mut))==1){
				top_mut_symbol<-sub("AMP_","",top_mut)
			}else if(length(grep("DEL_",top_mut))==1){
				top_mut_symbol<-sub("DEL_","",top_mut)
			}else if(length(grep("MET_",top_mut))==1){
				top_mut_symbol<-sub("MET_","",top_mut)
			}
			top_mut_symbol<-entrez2gene[top_mut_symbol]
			cat(top_mut_symbol,"is the top mut symbol\n")
			top_mut_colorsamples<-colAnnotation[top_mut_symbol,]
			top_mut_colorsamples<-top_mut_colorsamples[top_mut_colorsamples!="white"]
			top_mut_colorsamples[top_mut_colorsamples=="black"]<-1
			top_mut_colorsamples[top_mut_colorsamples=="red"]<-1
			top_mut_colorsamples[top_mut_colorsamples=="blue"]<--1
			top_mut_samples<-as.numeric(top_mut_colorsamples)
			names(top_mut_samples)<-names(top_mut_colorsamples)
			
			# Calculate gsea and get the leading edge samples
			est <- gsea2.es.int(vipermat[seed,],top_mut_samples,1)
			leg <- names(est$rlist)[1:which.max(est$es)]
			leg <- leg[leg %in% names(top_mut_samples)]
			effectSamples<-leg
		} else{ # Alternatively to ledge, for this gene, get only samples with activity above or below a certain threshold
			if(mut_assoc_here[best_muts[1,1],best_muts[1,2]]>=activityThreshold){
				effect<-"positive"
			} else {
				effect<-"negative"
			}
			if(effect=="positive"){
				effectSamples<-names(vipermat_here[seed,vipermat_here[seed,]>=activityThreshold])
			} else if (effect=="negative") {
				effectSamples<-names(vipermat_here[seed,vipermat_here[seed,]<=-(activityThreshold)])
			}
		}
		
		# The cluster is restricted to samples where the seed activity has the correct sign, and where mutations are present
		clusterSamples<-intersect(effectSamples,mutationSamples)
		#colAnnotation<-colAnnotation[,clusterSamples]

		# At this point, a failsafe measure is necessary, in the case clusterSamples is empty or having only 1 sample
		if(length(clusterSamples)<=1){
			mut_pvalue_here<-mut_pvalue_here[!rownames(mut_pvalue_here)%in%seed,!colnames(mut_pvalue_here)%in%names(these_muts)]
			mut_assoc_here<-mut_assoc_here[!rownames(mut_assoc_here)%in%seed,!colnames(mut_assoc_here)%in%names(these_muts)]
			skippedStep<-TRUE
			next
		}

		
#		# Clustering based on mutations
#		num<-matrix(0,nrow=nrow(colAnnotation),ncol=ncol(colAnnotation))
#		colnames(num)<-colnames(colAnnotation)
#		rownames(num)<-rownames(colAnnotation)
#		for(i in 1:nrow(num)){
#			for (j in 1:ncol(num)){
#				num[i,j]<-mut2number(colAnnotation[i,j])
#			}
#		}
#		#colDendro<-hclust(as.dist(1-cor(num,method="pearson")),method="average")
		
		
		# This group vs. rest of the samples VIPER t-test: obtaining the top MRs defining this cluster
		seedmat<-vipermat[,clusterSamples]
		notseedmat<-vipermat[,setdiff(colnames(vipermat),clusterSamples)]
		ttests<-f.rtt.na(seedmat,notseedmat)$p.value
		mrs_p<-ttests[,1]
		mrs_p<-sort(mrs_p)
		ps_sig<-length(mrs_p[mrs_p<=0.05])
		if(ps_sig>maxmrs){
			ps_sig<-maxmrs
		}
		mrs<-names(mrs_p[1:ps_sig])
		
		# Put our seed on top and order by average activity
		mrs<-setdiff(mrs,seed)
		plotMatrix<-vipermat[mrs,]
		plotMatrix<-plotMatrix[order(apply(plotMatrix[,clusterSamples],1,mean),decreasing=TRUE),]
		plotMatrix<-rbind(vipermat[seed,],plotMatrix)
		rownames(plotMatrix)[1]<-seed
		rownames(plotMatrix)<-entrez2gene[rownames(plotMatrix)]
		
		# Plot
		jpeg(paste(sep="",tumAcro,"/clusters/cluster_",step,"_",entrez2gene[seed],".jpeg"),width=4000,height=4000,pointsize=40)
		source("scripts/fedemap.R")
		fedemap(plotMatrix=plotMatrix,
				relevantColumns=c(clusterSamples), #,"GT010","GT011"),
				main=paste(sep="",tumAcro," subcluster ",step," - ",entrez2gene[seed]," (p=",signif(seedp,3),")"),
				colAnnotation=colAnnotation
		)
		dev.off()
		
		
		### Remove used mrs
		vipermat_here<-vipermat_here[!rownames(vipermat_here)%in%mrs,]
		
		### Remove used mutation events
		mut_pvalue_here<-mut_pvalue_here[!rownames(mut_pvalue_here)%in%mrs,!colnames(mut_pvalue_here)%in%names(these_muts)]
		mut_assoc_here<-mut_assoc_here[!rownames(mut_assoc_here)%in%mrs,!colnames(mut_assoc_here)%in%names(these_muts)]

		### Remvoe used samples (completely not overlapping samples)
		#unassigned_samples<-setdiff(colnames(vipermat_here),colnames(seedmat))
		#vipermat_here<-vipermat_here[!rownames(vipermat_here)%in%mrs,unassigned_samples]
		
		
		
		# Fill another element of the clusterlist
		integratedActivity<-apply(vipermat[,clusterSamples],1,stouffer)
		
		clusterlist[step][[1]]<-list(
				seed=seed,
				samples=clusterSamples,
				mutations=names(top_muts),
				pvalue=seedp,
				mrs=mrs,
				integratedActivity=integratedActivity
		)
	}
	
	
	
	return(clusterlist)
}

