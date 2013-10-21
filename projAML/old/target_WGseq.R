data.cli <- read.table("ftp://caftpd.nci.nih.gov/pub/dcc_target/AML/clinical/TARGET_AML_ClinicalDataSet_20120702.txt", header = TRUE, sep = "\t", na.strings="NA")
# "USI"                                       "Study"                                    
# "Gender"                                    "Race"                                     
# "Ethnicity"                                 "Age.in.days"                              
# "WBC..x10.3.MicroLiter."                    "Bone.marrow.leukemic.blast.percentage...."
# "Peripheral.blasts...."                     "Bone.Marrow.blast...in.relapse"           
# "Peripheral.blast...in.relapse"             "CNS.disease"                              
# "Chloroma"                                  "FAB"                                      
# "t.6.9."                                    "t.8.21."                                  
# "t.3.5..q25.q34."                           "t.6.11..q27.q23."                         
# "t.9.11..p22.q23."                          "t.10.11..p11.2.q23."                      
# "t.11.19..q23.p13.1."                       "inv.16."                                  
# "del5q"                                     "del7q"                                    
# "del9q"                                     "monosomy.5"                               
# "monosomy.7"                                "trisomy.8"                                
# "trisomy.21"                                "MLL"                                      
# "Other.MLL"                                 "Minus.Y"                                  
# "Cytogenetic.Code.Other"                    "Complexity"                               
# "ISCN"                                      "FLT3.ITD.positive."                       
# "FLT3.ITD.allelic.ratio"                    "NPM.mutation."                            
# "CEBPA.mutation."                           "WT1.mutation."                            
# "c.Kit.mutation."                           "c.Kit.Exon.mutation"                      
# "MRD.at.end.of.course.1"                    "MRD...at.end.of.course.1"                 
# "MRD.at.end.of.course.2"                    "MRD...at.end.of.course.2"                 
# "CR.status.at.end.of.course.1"              "CR.status.at.end.of.course.2"             
# "Risk.group"                                "EFS.time..days."                          
# "EFS.event.type.ID"                         "OS.event.ID"                              
# "OS.time..days."                            "SCT.in.1st.CR" 
data.cli$FAB.char <- gsub("[^[A-Za-z/d]]","",data.cli$FAB)
data.cli$Age.yrs.n <- ceiling(as.numeric(as.character(data.cli$Age.in.days))/365)
data.cli$Age.in.days.n <- as.numeric(as.character(data.cli$Age.in.days))
data.cli$OS.time..days.n <- as.numeric(as.character(data.cli$OS.time..days.))
class(data.cli$OS.time..days.)
hist(data.cli$Age.yrs.n)
data.cli$Bone.marrow.leukemic.blast.n<- as.numeric(as.character(data.cli$Bone.marrow.leukemic.blast.percentage....))
summary(data.cli$Bone.marrow.leukemic.blast.n)


install.packages("ggplot2")
require("ggplot2")

qplot(x=log(Age.in.days.n), y =log(OS.time..days.n), data= data.cli, colour=factor(FAB))
qplot(x=Age.yrs.n, y =OS.time..days.n, data= data.cli.osd1, colour=factor(FAB.char))+lines(x=14,col ="red")
qplot(y=OS.time..days.n,x=Age.in.days.n+Bone.marrow.leukemic.blast.n, data= data.cli, colour=factor(FAB))
qplot(y=OS.time..days.n,x=Bone.marrow.leukemic.blast.n, data= data.cli)


data.cli.osd1 <- data.cli[which(data.cli$OS.time..days.n<=2000),]
data.cli.osd1 <- data.cli.osd1[grep("[[:digit:]]",data.cli.osd1$FAB.char),]

qplot(y=log(OS.time..days.n),x=log(Bone.marrow.leukemic.blast.n), data= data.cli.osd1, colour= factor(FAB.char))



len <- dim(data.cli)[2]
for (i in 1: len) 
 grep("M0149",data.cli[,i])
