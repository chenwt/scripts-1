# TODO: Add comment
# 
# Author: fmg2117
###############################################################################
# This jewel converts a snv mutation matrix into a regulon object
convert2regulon<-function(x,sign=NULL) {
	x <- sign(tapply(x, names(x), mean))
	if(sign=="plus"){
		x <- x[x>0]
	} else if (sign=="minus"){
		x <- x[x<0]
	} else if (sign=="both") {
		x<-x[x!=0]
	}
	list(tfmode=x, likelihood=rep(1, length(x)))
}
# This jewel converts a fcnv mutation list into a regulon object
convertCNV2regulon<-function(x,sign=NULL) {
	if(sign=="plus"){
		x <- x[x>0]
	} else if (sign=="minus"){
		x <- x[x<0]
	} else if (sign=="both") {
		x<-x[x!=0]
	}
	list(tfmode=x, likelihood=rep(1, length(x)))
}
mut2number<-function(x){
	if(x=="black"){
		return(1)
	} else if (x=="white"){
		return(0) 
	} else if (x=="red") {
		return(1) 
	} else if (x=="blue") {
		return(-1)
	} else if (x=="green") {
		return(1)
	}
}
z2p<-function(z){
	pnorm(abs(z), lower.tail=F)*2
}
library(marina)

rootFolder<-"C:/Users/fmg2117/Dropbox/Projects_CUMC/panclustering"
#rootFolder<-"/ifs/scratch/c2b2/ac_lab/malvarez/panCancer"
tumAcro<-"stad"
setwd(rootFolder)
reportFile<-paste(sep="",tumAcro,"/",tumAcro,"-iterativeClustering_report.txt")

load("scripts/entrez2gene.rda")
source("scripts/panCancerFunctions.r")

### For a specific cancer, obtain:
# 1 - Functional SNVs
# 2 - Functional CNVs
# 3 - Functional Methylations

### Calculate activity to fMutation NESs

### Improvements
# 1 - Treat different signs in the activity/mutation matrix as different events
# 2 - Keep only mutations in the leading edge of the enrichment
# 3 - Robustness test: remove first gene and try again for consistency of the clustering; then remove first and second, etc.
# 4 - enrichment based on /data/malvarez/databases/kegg

cat(paste("Report for",tumAcro,"iterative clustering\n"),file=reportFile)

######################## Load data
# CCLE-TCGA normalized
load(paste("cclePCA/normexpressions_",tumAcro,".rda",sep=""))
load(paste("cclePCA/vipers_",tumAcro,".rda",sep=""))
load(paste("cclePCA/aggregatedVipers_",tumAcro,".rda",sep=""))
vipermat_aggregated<-vipermat
colnames(vipermat_aggregated)[grep("^TCGA",colnames(vipermat_aggregated),perl=TRUE)]<-substr(colnames(vipermat_aggregated)[grep("^TCGA",colnames(vipermat_aggregated),perl=TRUE)],1,12)
rm(vipermat)

# Standard Expression
if(!file.exists(paste(tumAcro,"/",tumAcro,"_withPatients_rnaseq-counts-normalized.rda",sep=""))){
	load(paste(tumAcro,"/",tumAcro,"_withPatients_rnaseq-counts.rda",sep=""))
	expmat<-DEtransform(dset)
	dim(expmat) # 22091 307
	rm(dset)
	save(expmat,samples,file=paste(tumAcro,"/",tumAcro,"_withPatients_rnaseq-counts-normalized.rda",sep=""))
} else {
	load(paste(tumAcro,"/",tumAcro,"_withPatients_rnaseq-counts-normalized.rda",sep=""))
}
normexpmat<-expmat[, as.character(samples[ which(samples[,2]=="Normal"), 1]) ]
tumorexpmat<-expmat[, as.character(samples[ which(samples[,2]=="Tumor"), 1]) ]
colnames(normexpmat)<-substr(colnames(normexpmat),1,12)
colnames(tumorexpmat)<-substr(colnames(tumorexpmat),1,12)


# Regulons
load(paste(tumAcro,"/",tumAcro,"_tcga_rnaseq_tf-regulon.rda",sep=""))
regul_tf<-regul
load(paste(tumAcro,"/",tumAcro,"_tcga_rnaseq_sig-regulon.rda",sep=""))
regul_sig<-regul
regul <- c(regul_tf, regul_sig[!(names(regul_sig) %in% names(regul_tf))])


### Standard VIPER Tumor vs. normal
if(!file.exists(paste(tumAcro,"/",tumAcro,"_withPatients-viper-results.rda",sep=""))){
	# T-statistics signature: TCGA tumor samples vs TCGA normal samples
	tss <- apply(tumorexpmat,2,
			function(x,normexpmat){
				tmp <- f.rtt.na(x - normexpmat)
				return((qnorm(tmp$p.value/2, lower.tail=F)*sign(tmp$statistic))[, 1])
			}, normexpmat=normexpmat)
	# Null model for Tumor vs Normal contrast
	set.seed(1453)
	tssnull <- tumorNormalNull(tumorexpmat, normexpmat)
	# Good old VIPER for Tumor vs. Normal TCGA signatures
	vipermat <- viper(tss,regul,tssnull,method="none", minsize=20)
	save(vipermat,file=paste(tumAcro,"/",tumAcro,"_withPatients-viper-results.rda",sep=""))
	
} else {
	load(paste(tumAcro,"/",tumAcro,"_withPatients-viper-results.rda",sep=""))
}


# Annotations
load(paste(sep="",tumAcro,"/",tumAcro,"_withPatients-snp.rda"))
load(paste(sep="",tumAcro,"/",tumAcro,"_withPatients-fcnv.rda"))
load(paste(sep="",tumAcro,"/",tumAcro,"-clinical.rda"))
load(paste(sep="",tumAcro,"/",tumAcro,"-fmeth.rda"))

### Formatting Mutation data
# Convert SNV data into a regulon object
if (!file.exists(paste(sep="",tumAcro,"/",tumAcro,"-snp_regulon.rda"))) {
	snp_regulon <- apply(rawsnp, 1, convert2regulon, sign="both")
	save(snp_regulon,file=paste(sep="",tumAcro,"/",tumAcro,"-snp_regulon.rda"))
} else {
	load(paste(sep="",tumAcro,"/",tumAcro,"-snp_regulon.rda"))
}

# Convert fCNV data into a regulon object
if (!file.exists(paste(sep="",tumAcro,"/",tumAcro,"-fcnv_regulon.rda"))) {
	fcnv_regulon_amp <- lapply(rawfcnv, convertCNV2regulon, sign="plus")
	fcnv_regulon_del <- lapply(rawfcnv, convertCNV2regulon, sign="minus")
	fcnv_regulon <- lapply(rawfcnv, convertCNV2regulon, sign="both")
	save(fcnv_regulon,fcnv_regulon_amp,fcnv_regulon_del, file=paste(sep="",tumAcro,"/",tumAcro,"-fcnv_regulon.rda"))
} else {
	load(paste(sep="",tumAcro,"/",tumAcro,"-fcnv_regulon.rda"))
}

# Convert methylation data into a regulon object
if (!file.exists(paste(sep="",tumAcro,"/",tumAcro,"-fmeth_regulon.rda"))) {
	fmeth_regulon <- lapply(rawfmeth, convertCNV2regulon,sign="both")
	save(fmeth_regulon,file=paste(sep="",tumAcro,"/",tumAcro,"-fmeth_regulon.rda"))
} else {
	load(paste(sep="",tumAcro,"/",tumAcro,"-fmeth_regulon.rda"))
}
pdf(paste(sep="",tumAcro,"/plots/meth_number_density.pdf"))
plot(density(unlist(lapply(fmeth_regulon,function(x){length(x$tfmode)}))))
dev.off()

# How many samples have mutations/methylations?
samples_with_snp<-unique(unlist(lapply(snp_regulon,function(x){names(x$tfmode)})))
samples_with_fcnv<-unique(unlist(lapply(fcnv_regulon,function(x){names(x$tfmode)})))
samples_with_fmeth<-unique(unlist(lapply(fmeth_regulon,function(x){names(x$tfmode)})))
cat(
		paste(sep="",
				"Total normal samples in ",tumAcro,": ",ncol(normexpmat),"\n",
				"Total tumor samples in ",tumAcro,": ",ncol(vipermat),"\n",
				"This number of ",tumAcro," samples have reported SNVs: ",length(samples_with_snp),"\n",
				"This number of ",tumAcro," samples have reported fCNVs: ",length(samples_with_fcnv),"\n",
				"This number of ",tumAcro," samples have functional methylations above beta threshold ",meth_threshold,": ",length(samples_with_fmeth),"\n",
				"This number of ",tumAcro," samples have reported mutations: ",length(union(samples_with_snp,samples_with_fcnv)),"\n",
				"This number of ",tumAcro," samples have reported mutations/methylations: ",length(union(samples_with_snp,union(samples_with_fcnv,samples_with_fmeth))),"\n\n"
				),
		file=reportFile,
		append=TRUE
)


###################### ENRICHMENT ANALYSIS
### SNV Enrichment analysis
minsamples<-10
if (!file.exists(paste(sep="",tumAcro,"/",tumAcro,"-snp-association.rda"))) {
	snp_assoc <- viper(t(vipermat), snp_regulon, tw=0, lw=0, method="none", eset.filter=F, minsize=minsamples) # 1714 genes with >10 SNVs
	save(snp_assoc, file=paste(sep="",tumAcro,"/stad-snp-association.rda"))
} else {
	load(paste(sep="",tumAcro,"/",tumAcro,"-snp-association.rda"))
}

### CNV Enrichment analysis
minsamples<-10
if (!file.exists(paste(sep="",tumAcro,"/",tumAcro,"-fcnv-association.rda"))) {
	fcnv_assoc_amp <- viper(t(vipermat), fcnv_regulon_amp, tw=0, lw=0, method="none", eset.filter=F, minsize=minsamples) 
	fcnv_assoc_del <- viper(t(vipermat), fcnv_regulon_del, tw=0, lw=0, method="none", eset.filter=F, minsize=minsamples) 
	fcnv_assoc <- viper(t(vipermat), fcnv_regulon, tw=0, lw=0, method="none", eset.filter=F, minsize=minsamples) 
	save(fcnv_assoc,fcnv_assoc_amp,fcnv_assoc_del, file=paste(sep="",tumAcro,"/",tumAcro,"-fcnv-association.rda"))
} else {
	load(paste(sep="",tumAcro,"/",tumAcro,"-fcnv-association.rda"))
}

### Meth Enrichment analysis (SWITCH THE SIGN!)
minsamples<-10
if (!file.exists(paste(sep="",tumAcro,"/",tumAcro,"-fmeth-association.rda"))) {
	fmeth_assoc <- viper(-t(vipermat), fmeth_regulon, tw=0, lw=0, method="none", eset.filter=F, minsize=minsamples) 
	save(fmeth_assoc, file=paste(sep="",tumAcro,"/",tumAcro,"-fmeth-association.rda"))
} else {
	load(paste(sep="",tumAcro,"/",tumAcro,"-fmeth-association.rda"))
}


### Format the association tables
snp_assoc<-t(snp_assoc)
fcnv_assoc<-t(fcnv_assoc)
fcnv_assoc_amp<-t(fcnv_assoc_amp)
fcnv_assoc_del<-t(fcnv_assoc_del)
fmeth_assoc<-t(fmeth_assoc)

### Convert NESs into p-values
snp_pvalue<-matrix(p.adjust(z2p(snp_assoc),method="fdr"),ncol=ncol(snp_assoc),nrow=nrow(snp_assoc))
rownames(snp_pvalue)<-rownames(snp_assoc)
colnames(snp_pvalue)<-colnames(snp_assoc)

fcnv_pvalue<-matrix(p.adjust(z2p(fcnv_assoc),method="fdr"),ncol=ncol(fcnv_assoc),nrow=nrow(fcnv_assoc))
rownames(fcnv_pvalue)<-rownames(fcnv_assoc)
colnames(fcnv_pvalue)<-colnames(fcnv_assoc)

fcnv_pvalue_amp<-matrix(p.adjust(z2p(fcnv_assoc_amp),method="fdr"),ncol=ncol(fcnv_assoc_amp),nrow=nrow(fcnv_assoc_amp))
rownames(fcnv_pvalue_amp)<-rownames(fcnv_assoc_amp)
colnames(fcnv_pvalue_amp)<-colnames(fcnv_assoc_amp)

fcnv_pvalue_del<-matrix(p.adjust(z2p(fcnv_assoc_del),method="fdr"),ncol=ncol(fcnv_assoc_del),nrow=nrow(fcnv_assoc_del))
rownames(fcnv_pvalue_del)<-rownames(fcnv_assoc_del)
colnames(fcnv_pvalue_del)<-colnames(fcnv_assoc_del)

fmeth_pvalue<-matrix(p.adjust(z2p(fmeth_assoc),method="fdr"),ncol=ncol(fmeth_assoc),nrow=nrow(fmeth_assoc))
rownames(fmeth_pvalue)<-rownames(fmeth_assoc)
colnames(fmeth_pvalue)<-colnames(fmeth_assoc)


### Combine together the mutational (methylation if available) data
mut_assoc<-cbind(snp_assoc,fcnv_assoc_amp,fcnv_assoc_del,fmeth_assoc)
colnames(mut_assoc)<-c(
		paste("SNV_",colnames(snp_assoc),sep=""),
		paste("AMP_",colnames(fcnv_assoc_amp),sep=""),
		paste("DEL_",colnames(fcnv_assoc_del),sep=""),
		paste("MET_",colnames(fmeth_assoc),sep="")
)
### P-value conversion, again (BH comparable now)
mut_pvalue<-matrix(p.adjust(pnorm(abs(mut_assoc), lower.tail=F)*2,method="fdr"),ncol=ncol(mut_assoc),nrow=nrow(mut_assoc))
rownames(mut_pvalue)<-rownames(mut_assoc)
colnames(mut_pvalue)<-colnames(mut_assoc)
# Significance plot
jpeg(paste(tumAcro,"/plots/",tumAcro,"-significanceDistribution.jpeg",sep=""),width=1500,height=1500,pointsize=30)
plot(density(mut_pvalue[,grep("MET_",colnames(mut_pvalue))]),
		lwd=2,col="green",
		main="Significance of association between VIPER activity and mutations",
		xlab="P-value"
)
mtext("Stomach Adenocarcinoma - FDR correction")
abline(v=0.05,col="black",lty=2)
lines(density(mut_pvalue[,grep("AMP_",colnames(mut_pvalue))]),
		lwd=2,col="red",
)
lines(density(mut_pvalue[,grep("SNV_",colnames(mut_pvalue))]),
		lwd=2,col="black",
)
lines(density(mut_pvalue[,grep("DEL_",colnames(mut_pvalue))]),
		lwd=2,col="blue",
)
legend("top",legend=c("SNV","CNV amp","CNV del","MET"),col=c("black","red","blue","green"),lty=1,lwd=2)
dev.off()



##########################################################
################### ITERATIVE CLUSTERING
source(paste(sep="","scripts/loopClustering.R"))
clusterlist<-loopClustering(
		vipermat=vipermat,
		mut_pvalue=mut_pvalue,
		mut_assoc=mut_assoc,
		entrez2gene=entrez2gene,
		tumorexpmat=tumorexpmat,
		normexpmat=normexpmat,
		regul=regul,
		maxstep=50,
		maxmrs=20,
		maxmutations=20,
		ledge=FALSE,
		activityThreshold=0.84 # 0.84 means 20% chance of being an increased/decreased activity by chance over normal
)
assigned_patients<-unique(unlist(lapply(clusterlist,function(x){x$samples})))
unassigned_patients<-setdiff(colnames(vipermat),assigned_patients)
cat(
		paste(sep="",
				"Total normal samples in ",tumAcro,": ",ncol(normexpmat),"\n"
		),
		file=reportFile,
		append=TRUE
)


### TODO TODO wrong wrong assignment of samples with DELETIONS!!! (e.g. cluster 11: OR5AN1)

### Show the number of samples added after each step
# and the number of total samples
# and the number of samples with >10 mutations
samples_all<-colnames(vipermat)
samples_with_snp
samples_with_fcnv
samples_with_mut<-union(samples_with_snp,samples_with_fcnv)

samples_included<-c()
mrs_included<-c()
nr_samples_included<-c()
nr_mrs_included<-c()
length_clusters<-c()
for(i in 1:length(clusterlist)){
	newsamples<-clusterlist[i][[1]]$samples
	newmrs<-clusterlist[i][[1]]$mrs
	samples_included<-union(samples_included,newsamples)
	mrs_included<-union(mrs_included,newmrs)
	nr_samples_included<-c(nr_samples_included,length(samples_included))
	nr_mrs_included<-c(nr_mrs_included,length(newmrs))
	length_clusters<-c(length_clusters,length(clusterlist[i][[1]]$samples))
}




pdf(paste(sep="",tumAcro,"/plots/lineplot_includedSamples.pdf"),width=10,height=6)
plot(nr_samples_included,lwd=2,type="l",
		xlab="clustering step",ylab="Number of samples (or MRs)",
		main=paste(sep="",tumAcro," - included samples during first clustering iteration"),
		ylim=c(0,length(samples_all)),
		col=1

)
lines(length_clusters,lwd=2,lty=2,col=2)
lines(nr_mrs_included,lwd=2,lty=2,col=3)

abline(h=length(samples_with_mut),lty=2,col="darkgrey")
text(x=1,y=length(samples_with_mut),labels="Nr. samples with reported mutations",pos=4)
abline(h=length(samples_all),lty=2,col="darkgrey")
text(x=1,y=length(samples_all),labels="Nr. samples with reported expression",pos=4)

legend("topright",
		legend=c("Total samples included","Size of new cluster introduced","New MRs introduced in each step"),
		col=c(1:3),lty=c(1,2,2),lwd=2,
		bg="white"
)
dev.off()

# P-value trend
ps<-c()
for(i in 1:length(clusterlist)){
	ps<-c(ps,clusterlist[i][[1]]$pvalue)
}
pdf(paste(sep="",tumAcro,"/plots/lineplot_pvalues.pdf"),width=10,height=6)
plot(ps,xlab="step",ylab="P-value of activity/mutation used as cluster seed",
		main=paste(sep="",tumAcro," - p-values of iterative clustering steps")
		)
dev.off()




### Clustering: are the unassigned patients making other clusters? 
library(tsne)
ttt<-tsne(as.dist(signatureDistance(vipermat)))
rownames(ttt)<-colnames(vipermat)

pdf(paste(sep="",tumAcro,"/plots/tsne_afterIterativeClustering.pdf"),width=8,height=8)
plot(ttt,col="white",main=paste(tumAcro,"- tsne on VIPER vs. Normal | Signature Distance"))
points(
		ttt[assigned_patients,1],
		ttt[assigned_patients,2],
		pch=1)
points(
		ttt[grep("GT010",rownames(ttt)),1],
		ttt[grep("GT010",rownames(ttt)),2],
		pch=2)
points(
		ttt[grep("GT011",rownames(ttt)),1],
		ttt[grep("GT011",rownames(ttt)),2],
		pch=3)
points(
		ttt[unassigned_patients,1],
		ttt[unassigned_patients,2],
		pch=4)
legend("topleft",
		col=c("black"),
		pch=c(1,2,3,4),
		legend=c("Assigned","Cornell","Columbia","Not assigned")
)
dev.off()


### Intersection: overlap between clusters
# Top 5 venn
library(VennDiagram)
venngroups<-lapply(clusterlist,function(x){x$samples})
names(venngroups)<-c(1:length(venngroups))
venn.diagram(venngroups[1:5],fill = c(1:5),
		alpha = c(0.5), cex = 2,cat.fontface = 4,lty =2, 
		filename = paste(sep="",tumAcro,"/plots/venn_top5.tif"));

# Pairwise overlap between groups
library(plotrix)
sample_overlap_integer<-matrix(0,ncol=length(clusterlist),nrow=length(clusterlist))
for(i in 1:length(clusterlist)){
	for(j in 1:length(clusterlist)){
		o<-length(intersect(clusterlist[i][[1]]$samples,clusterlist[j][[1]]$samples))
		sample_overlap_integer[i,j]<-o
	}
}
sample_overlap_relative<-matrix(0.0,ncol=length(clusterlist),nrow=length(clusterlist))
for(i in 1:length(clusterlist)){
	for(j in 1:length(clusterlist)){
		o<-length(intersect(clusterlist[i][[1]]$samples,clusterlist[j][[1]]$samples))/length(clusterlist[i][[1]]$samples)
		sample_overlap_relative[i,j]<-o
	}
}
#rownames(sample_overlap_relative)<-entrez2gene[unlist(lapply(clusterlist,function(x){x$seed}))]
pdf(paste(sep="",tumAcro,"/plots/sample_overlap.pdf"),width=12,height=18)
par(mfrow=c(2,1))
color2D.matplot(sample_overlap_integer,show.values=TRUE,main="Cluster overlap as number of samples",vcex=0.8)
color2D.matplot(sample_overlap_relative,show.values=TRUE,main="Cluster overlap, relative to row cluster size",vcex=0.8)
dev.off()



##############################################################################################
##############################################################################################
##############################################################################################
###### Robustness: remove the best cluster and check if the clustering remains the same ######
# Calculate global pairwise signatureDistance, PCC and SCC
dist_signdist<-signatureDistance(vipermat,ws=2)
dist_pearson<-cor(vipermat)
dist_spearman<-cor(vipermat,method="s")

### 1) remove one by one the mutated samples and see where they fall (signature distance distribution)
### 2) best discriminating distance metrics (pearson, spearman, signature distance with weight 1, 2 or 3)
# Rank distribution for robust cluster assignation
robustness_ranks_signdist<-list()
robustness_ranks_spearman<-list()
for(i_cluster in 1:length(clusterlist)){
	this_cluster<-clusterlist[i_cluster][[1]]
	for (j_sample in 1:length(this_cluster$samples)){
		this_sample<-this_cluster$samples[j_sample]
		# For each sample, calculate the distance vs. all the other samples
		this_signdist<-list()
		this_spearman<-list()
		for(k_allcluster in 1:length(clusterlist)){
			this_signdist[k_allcluster][[1]]<-dist_signdist[setdiff(clusterlist[k_allcluster][[1]]$samples,this_sample),this_sample]
			this_spearman[k_allcluster][[1]]<-dist_spearman[setdiff(clusterlist[k_allcluster][[1]]$samples,this_sample),this_sample]
		}
		# Average the distance # BEWARE: signature Distance is actually signature SIMILARITY (e.g. positive is similar)
		avg_signdist<-unlist(lapply(this_signdist,mean))
		rank_signdist<-rank(avg_signdist)
		avg_spearman<-unlist(lapply(this_spearman,mean))
		rank_spearman<-rank(avg_spearman)
		# Add the rank of the same cluster to the cluster robustness
		robustness_ranks_signdist[i_cluster][[1]]<-c(robustness_ranks_signdist[i_cluster][[1]] , rank_signdist[i_cluster])
		robustness_ranks_spearman[i_cluster][[1]]<-c(robustness_ranks_spearman[i_cluster][[1]] , rank_spearman[i_cluster])
	}
}

# Random rank distribution
random_robustness_ranks_signdist<-list()
for(i in 1:1){
	for(i_cluster in 1:length(clusterlist)){
		this_cluster<-clusterlist[i_cluster][[1]]
		for (j_sample in 1:length(this_cluster$samples)){
			this_sample<-sample(assigned_patients,1) # Here's the randomization
			# For each sample, calculate the distance vs. all the other samples
			this_signdist<-list()
			for(k_allcluster in 1:length(clusterlist)){
				this_signdist[k_allcluster][[1]]<-dist_signdist[setdiff(clusterlist[k_allcluster][[1]]$samples,this_sample),this_sample]
			}
			# Average the distance # BEWARE: signature Distance is actually signature SIMILARITY (e.g. positive is similar)
			avg_signdist<-unlist(lapply(this_signdist,mean))
			rank_signdist<-rank(avg_signdist)
			# Add the rank of the same cluster to the cluster robustness
			random_robustness_ranks_signdist[i_cluster][[1]]<-c(random_robustness_ranks_signdist[i_cluster][[1]] , rank_signdist[i_cluster])
		}
	}
}


pdf(paste(sep="",tumAcro,"/plots/boxplots_robustness_ranks_signdist.pdf"),width=12,height=6)
boxplot(robustness_ranks_signdist,ylab="Rank",main=paste(sep="",length(clusterlist)," clusters in ",tumAcro),cex.axis=0.6)
mtext("Signature distance, ws=2")
dev.off()

pdf(paste(sep="",tumAcro,"/plots/boxplots_robustness_ranks_spearman.pdf"),width=12,height=6)
boxplot(robustness_ranks_spearman,ylab="Rank",main=paste(sep="",length(clusterlist)," clusters in ",tumAcro),cex.axis=0.6)
mtext("Spearman correlation")
dev.off()

pdf(paste(sep="",tumAcro,"/plots/boxplots_random_robustness_ranks_signdist.pdf"),width=12,height=6)
boxplot(random_robustness_ranks_signdist,ylab="Rank",main=paste(sep="",length(clusterlist)," clusters in ",tumAcro),cex.axis=0.6)
mtext("Signature distance, ws=2")
dev.off()

pdf(paste(sep="",tumAcro,"/plots/density_robustness_ranks.pdf"),width=6,height=6)
plot(density(unlist(robustness_ranks_signdist)),lwd=2,col="black",lty=2)
lines(density(unlist(random_robustness_ranks_signdist)),lwd=2,col="red",lty=1)
lines(density(unlist(robustness_ranks_spearman)),lwd=3,col="blue",lty=3)
legend("topleft",legend=c("Signature Distance rank","Randomized SD rank","Spearman Rank"),lty=c(2,1,3),col=c("black","red","blue"),lwd=2)
dev.off()

pdf(paste(sep="",tumAcro,"/plots/boxplot_robustness_ranks.pdf"),width=6,height=6)
boxplot(unlist(random_robustness_ranks_signdist),
		unlist(robustness_ranks_signdist),
		unlist(robustness_ranks_spearman),
		names=c("Random SD ws=2","SD ws=2","Spearman"),
		main="Rank assignment of sample on original cluster"
)
mtext("Higher is better")
dev.off()




### TODO: 3) remove the starting point and check the consistency of generated clusters



#########################################################################################
### Assignment: based on the signatureDistance matrix, ASSIGN THE UNASSIGNED PATIENTS ###
#########################################################################################
clusterlist2<-clusterlist
nr_newlyassigned<-rep(0,length(clusterlist))
for (this_sample in unassigned_patients){
	# For each sample, calculate the distance vs. all the clusters
	this_signdist<-list()
	for(k_allcluster in 1:length(clusterlist)){
		this_signdist[k_allcluster][[1]]<-dist_signdist[clusterlist[k_allcluster][[1]]$samples,this_sample]
	}
	# Average the distance # BEWARE: signature Distance is actually signature SIMILARITY (e.g. positive is similar)
	avg_signdist<-unlist(lapply(this_signdist,mean))
	# Assign the sample to the top ranking cluster
	matchcluster<-which.max(avg_signdist)
	clusterlist2[matchcluster][[1]]$newsamples<-c(clusterlist2[matchcluster][[1]]$newsamples,this_sample)
	nr_newlyassigned[matchcluster]<-nr_newlyassigned[matchcluster]+1
}
pdf(paste(sep="",tumAcro,"/plots/barplot_nr_newlyassigned.pdf"),width=6,height=6)
barplot(nr_newlyassigned,names.arg=1:length(nr_newlyassigned),main="Assignment of unassigned samples per cluster",ylab="Number of new samples",xlab="Cluster",cex.names=0.6)
dev.off()

### Cluster Plots with newly assigned samples:
for(i in 1:length(clusterlist2)) {
	step<-i
	seed<-clusterlist2[i][[1]]$seed
	seedp<-clusterlist2[i][[1]]$pvalue
	mrs<-clusterlist2[i][[1]]$mrs
	top_muts<-clusterlist2[i][[1]]$mutations
	
	oldClusterSamples<-clusterlist2[i][[1]]$samples
	newClusterSamples<-clusterlist2[i][[1]]$newsamples
	clusterSamples<-c(oldClusterSamples,newClusterSamples)
	
	# Prepare the plot matrix
	plotMatrix<-vipermat[mrs,]
	plotMatrix<-plotMatrix[order(apply(plotMatrix[,clusterSamples],1,mean),decreasing=TRUE),]
	plotMatrix<-rbind(vipermat[seed,],plotMatrix)
	rownames(plotMatrix)[1]<-seed
	rownames(plotMatrix)<-entrez2gene[rownames(plotMatrix)]
	
	# Generate column association (mutation) color object
	colAnnotation<-matrix("white",nrow=length(top_muts),ncol=ncol(vipermat))
	rownames(colAnnotation)<-top_muts
	colnames(colAnnotation)<-colnames(vipermat)
	for (mutgene in top_muts) {
		origene<-mutgene
		if(length(grep("SNV_",mutgene))==1){
			mutgene<-sub("SNV_","",mutgene)
			col<-rep("white",length(rawsnp[mutgene,]))
			col[rawsnp[mutgene,]!=0]<-"black"
			names(col)<-names(rawsnp[mutgene,])
		} else if(length(grep("AMP_",mutgene))==1){
			mutgene<-sub("AMP_","",mutgene)
			col<-rep("white",length(rawfcnv[[mutgene]]))
			col[rawfcnv[[mutgene]]>0]<-"red"
			col[rawfcnv[[mutgene]]<0]<-"blue"
			names(col)<-names(rawfcnv[[mutgene]])
		} else if(length(grep("DEL_",mutgene))==1){
			mutgene<-sub("DEL_","",mutgene)
			col<-rep("white",length(rawfcnv[[mutgene]]))
			col[rawfcnv[[mutgene]]>0]<-"red"
			col[rawfcnv[[mutgene]]<0]<-"blue"
			names(col)<-names(rawfcnv[[mutgene]])
		} else if(length(grep("MET_",mutgene))==1){
			mutgene<-sub("MET_","",mutgene)
			col<-rep("white",length(rawfmeth[[mutgene]]))
			col[rawfmeth[[mutgene]]!=0]<-"green"
			names(col)<-names(rawfmeth[[mutgene]])
		}
		commonsamples<-intersect(names(col),colnames(colAnnotation))
		colAnnotation[origene,commonsamples]<-col[commonsamples]
	}
	rownames(colAnnotation)<-paste(sep="",substr(rownames(colAnnotation),1,4),entrez2gene[substr(rownames(colAnnotation),5,nchar(rownames(colAnnotation)))])
	
	
	jpeg(paste(sep="",tumAcro,"/clusters/newcluster_",step,"_",entrez2gene[seed],".jpeg"),width=4000,height=2000,pointsize=40)
	source("scripts/fedemap.R")
	fedemap(plotMatrix=plotMatrix,
			relevantColumns=oldClusterSamples,
			main=paste(sep="",tumAcro," subcluster ",step," - ",entrez2gene[seed]," (p=",signif(seedp,3),")"),
			colAnnotation=colAnnotation,
			blueColumns=newClusterSamples
	)
	dev.off()
}


# Where is our patient?
for (i in 1:length(clusterlist2)){
	cat(i,"GT011"%in%clusterlist2[i][[1]]$newsamples,"\n")
}


#######################################################################################
################################ ASSIGN THE CELL LINES ################################
#######################################################################################
signdist_aggregated<-signatureDistance(vipermat_aggregated,ws=2)

## Test removal of inter-dataset variance:
#ttt<-tsne(as.dist(signdist_aggregated))
#rownames(ttt)<-colnames(vipermat_aggregated)
#plot(ttt,col="white",main=paste(tumAcro,"- tsne on VIPER vs. Normal | PCC"))
#points(
#		ttt[grep("TCGA",rownames(ttt)),1],
#		ttt[grep("TCGA",rownames(ttt)),2],
#		pch=1)
#points(
#		ttt[grep("GT010",rownames(ttt)),1],
#		ttt[grep("GT010",rownames(ttt)),2],
#		pch=2)
#points(
#		ttt[grep("_STOMACH",rownames(ttt)),1],
#		ttt[grep("_STOMACH",rownames(ttt)),2],
#		pch=4)
#legend("topleft",
#		col=c("black"),
#		pch=c(1,2,4),
#		legend=c("TCGA","Cornell","Cell Line")
#)


### Map cell lines to their best cluster
clusterlist3<-clusterlist2
nr_newlyassigned<-rep(0,length(clusterlist))
for (this_sample in colnames(viper_ccle)){
	# For each sample, calculate the distance vs. all the clusters
	this_signdist<-list()
	for(k_allcluster in 1:length(clusterlist)){
		patients_here<-c(clusterlist[k_allcluster][[1]]$samples,clusterlist[k_allcluster][[1]]$newsamples)
		this_signdist[k_allcluster][[1]]<-signdist_aggregated[patients_here,this_sample]
	}
	# Average the distance # BEWARE: signature Distance is actually signature SIMILARITY (e.g. positive is similar)
	avg_signdist<-unlist(lapply(this_signdist,mean))
	# Assign the sample to the top ranking cluster
	matchcluster<-which.max(avg_signdist)
	clusterlist3[matchcluster][[1]]$cclesamples<-c(clusterlist3[matchcluster][[1]]$cclesamples,this_sample)
	nr_newlyassigned[matchcluster]<-nr_newlyassigned[matchcluster]+1
}
pdf(paste(sep="",tumAcro,"/plots/barplot_nr_newlyassigned_ccle.pdf"),width=6,height=6)
barplot(nr_newlyassigned,names.arg=1:length(nr_newlyassigned),main="Assignment of unassigned cell lines per cluster",ylab="Number of new samples",xlab="Cluster",cex.names=0.6)
dev.off()


save(clusterlist,clusterlist2,clusterlist3,file=paste(sep="",tumAcro,"/",tumAcro,"-panclusterlists.rda"))





### GO enrichment of MRs of each cluster


### TODO: rerun with sample exclusion
### TODO: highlight best mutation matched
### TODO: fix CNV bug

