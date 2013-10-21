
MAplot <- function(tumor.file,normal.file,t.num,n.num,del_reg.file){
	require(ggplot2)
	###------------------------------ load data
	tumor <- read.table(tumor.file,sep="\t",header=T,stringsAsFactors=F)
	normal <- read.table(normal.file,sep="\t",header=T,stringsAsFactors=F)
	tumor$FPKM[tumor$FPKM == 0] <- 1e-25
	normal$FPKM[normal$FPKM == 0] <- 1e-25
	M <- log(tumor$FPKM,2)-log(normal$FPKM,2)
	A <- 1/2 * (log(tumor$FPKM,2) + log(normal$FPKM,2))
	data.plot <- data.frame(M=M,A=A)
	###------------------------------load deletion region data
	del.reg <- read.table(del_reg.file,sep="\t",header=T)
	del.reg <- del.reg[,2:4]
	colnames(del.reg) <- c("chromosome","start","end")
	del.reg$chromosome <- as.numeric(del.reg$chromosome)
	tmpvar <-  unlist(lapply(tumor$locus,FUN = function(x){unlist(strsplit(x,"-"))[1]}))
	emph.chr <-  unlist(lapply(tmpvar,FUN = function(x){unlist(strsplit(x,":"))[1]})) 
	emph.start <- as.numeric( unlist(lapply(tmpvar,FUN = function(x){unlist(strsplit(x,":"))[2]})) )
	emph.end <- as.numeric( unlist(lapply(tumor$locus,FUN = function(x){unlist(strsplit(x,"-"))[2]})))
	emph <- cbind(data.plot, emph.chr,emph.start,emph.end)
	colnames(emph) <- c(colnames(data.plot),"chr","start","end")

	##------------------------------
	emph.plot <- data.frame(M=numeric(),A=numeric(),chr=numeric(),start=numeric(),end=numeric())
	emph$chr <- as.numeric(emph$chr)

	for (i in 1:nrow(emph)) {
		for (j in 1:nrow(del.reg)){
			if( !is.na(emph[i,3]) & emph[i,3] == del.reg[j,1]){
				if( emph[i,4]>= del.reg[j,2] & emph[i,5] <= del.reg[j,3]){
			        emph.plot <- rbind(emph.plot,emph[i,])
			    	break
				}
			}
		next
		}
	}

	#------------------------------statistic test
	###------------------------------ plot 
	# colnames(data.plot) <- c("M","A")
	title <- paste("MAplot_", t.num,"_",n.num,sep="")
	pdf(paste(title,".pdf",sep=""))
	plot <- ggplot(aes(x=A, y=M), data=data.plot) + geom_point() + ylim(-100,100) + xlim(-60,15) 
	plot <- plot + geom_point(data=emph.plot,aes(x=A,y=M),color="red")
	plot <- plot + opts(panel.background = theme_rect(fill='white', color='black'))
	plot <- plot +  opts(title=title)
	print(plot)
	dev.off()

}







