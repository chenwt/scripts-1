#!/bin/Rscirpts
#Author: Jing He
#Date: Apr. 7th, 2013
#Last Update: 
#parameters: 
#example: Rscript 

###------------------------------header------------------------------
args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters")
	exit
}else{
	print(args)
}

###------------------------------start coding##------------------------------

if(Sys.info()["sysname"] == "Darwin") {
	setwd("/Volumes/ys_lab_scratch/jh3283/school/compGenomic/LAML/")
	print("work on IMAC!")
}else {
	setwd("~/SCRATCH/school/compGenomic/LAML/")
	print("Work On TITAN!")
}


##------------------------------load data
data <- read.delim("WGS/LAML_WGS.txt",header=T)
# [1] "analysis_id"             "state"                   "reason"                 
#  [4] "last_modified"           "upload_date"             "published_date"         
#  [7] "center_name"             "study"                   "aliquot_id"             
# [10] "sample_accession"        "legacy_sample_id"        "disease_abbr"           
# [13] "tss_id"                  "participant_id"          "sample_id"              
# [16] "analyte_code"            "sample_type"             "library_strategy"       
# [19] "platform"                "refassem_short_name"     "analysis_submission_uri"
# [22] "analysis_full_uri"       "analysis_data_uri"       "filesize"               
# [25] "filename"                "checksum"               
data2 <- read.delim("LAML_exome.txt",header=T)
dplot <- rbind(data[,c("filename","tss_id", "sample_type","disease_abbr","library_strategy","filesize")],
			data2[,c("filename","tss_id",  "sample_type","disease_abbr","library_strategy","filesize")])
dplot$sample_type <- as.character(dplot$sample_type)
dplot.m <- melt(dplot[,c("filename","sample_type","library_strategy","filesize")],
				id.var = "filename", 
				variable.name=c("sample_type"),
				measure.vars="filesize")
write.table(dplot,"LAML_WGS_Exome.txt",quote=NULL,sep="\t")
head(dplot.m)
attach(dplot)
plot(dplot[which(sample_type =="WGS"),"filesize"]/1e+09)




require(ggplot2)
g <- ggplot(dplot, aes(x=filename, y= filesize/1e+09, group=library_strategy)) + geom_line()
	# + opts(title = expression("TCGA LAML sequencing data summary"))

	# + labs(x = NULL, y = "File Size")
g
g <- g +  facet_grid(. ~ library_strategy)
g


# attach(dplot.m)
a <- ggplot(dplot, aes(x =as.character(sample_type), y = filesize/1e+09)) 
	# + opts(title = "TCGA LAML sequencing data summary") +
   # labs(x = NULL, y = "File Size",
   #     fill = NULL)
a + geom_histogram()

c <- b + facet_grid( ~ library_strategy.) + opts(legend.position = "none")
c



df.m <- melt(df)
df.m <- rename(df.m, c(X1 = "Period", X2 = "Region"))


