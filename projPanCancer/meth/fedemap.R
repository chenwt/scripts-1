#' FEDEMAP
#'
#' Heatmap with custom layout
#'
#' @param plotMatrix Numeric matrix, plotted as a blue/red heatmap, with defined column and row names
#' @param colDendro Hierarchical clustering object for columns
#' @param rowDendro Hierarchical clustering object for rows
#' @param colAnnotation matrix of colors (defined as color names or rgb) to plot on top of the heatmap
#' @param rowAnnotation
#' @param main The title of your plot
#' @param nbreaks Number of breaks in the red/blue heatmap (default: 30)
#' @return A plot
#' \describe{
#' \item{plot}{Plots several objects}
#' }
#' @export
require(gplots)

fedemap<-function(
		plotMatrix=NULL,
		blueColumns=NULL,
		relevantColumns=NULL,
		relevantRows=NULL,
		colDendro=NULL,
		rowDendro=NULL,
		colAnnotation=NULL,
		rowAnnotation=NULL,
		main="Title",
		nbreaks=40
){
	### Initialize parameters
	extreme=ceiling(max(abs(plotMatrix)))
	breaks <- seq(-extreme, extreme, length = nbreaks)
	ncol <- length(breaks) - 1
	col <- colorpanel(ncol,"blue","white","red")
	toprowOrder<-rev(order(plotMatrix[1,]))
	plotMatrix<-plotMatrix[,toprowOrder]
	
	
	### Define Layout
	layoutMatrix<-rbind(
			c(1,1,1),
			c(2,3,4),
			c(5,6,4),
			c(7,8,9)
	)
	
	layout(layoutMatrix,
			widths=c(1,15,1),
			heights=c(1,2,7,7)
	)
	#layout.show(n=12)
	par(mar=c(0,0,0,0))
	
	### 1 - Title
	plot(0,col="white",xlim=c(0,10),ylim=c(0,10),xaxt="n",yaxt="n",type="n",frame.plot=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
	text(5,5,labels=main,cex=2)
	
	### 2 - Key
	par(mar=c(1,2,1,0))
	image(t(as.matrix(breaks)),col=col,axes=FALSE,main="Key",cex.main=1,cex=1)
	mtext(text=pretty(breaks), side=2, line=0.3, at=(rank(pretty(breaks))-0.45)/length(pretty(breaks)), las=1,cex=0.6)
	par(mar=c(0,0,0,0))
	
	### 3 - Sample names
	plot(0,col="white",xlim=c(0,ncol(plotMatrix)),ylim=c(0,10),xaxt="n",yaxt="n",type="n",frame.plot=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
	text(1:ncol(plotMatrix),3,colnames(plotMatrix),pos=3,srt=90,col=
					ifelse(colnames(plotMatrix)%in%relevantColumns,"black","#00000000") # Transparency
	)
	text(1:ncol(plotMatrix),3,colnames(plotMatrix),pos=3,srt=90,
			col=ifelse(colnames(plotMatrix)%in%blueColumns,"blue","#00000000") # Transparency
	)
	
	
	### 4 - Empty
	plot(0,col="white",xlim=c(0,10),ylim=c(0,10),xaxt="n",yaxt="n",type="n",frame.plot=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
	
	### 5 - Mutation names
	par(mar=c(1,0,0,0))
	plot(0,col="white",xlim=c(0,10),ylim=c(0,nrow(colAnnotation)),xaxt="n",yaxt="n",type="n",frame.plot=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
	text(rep(10),(1:nrow(colAnnotation)-0.5),labels=rev(rownames(colAnnotation)),cex=1,pos=2)
	par(mar=c(0,0,0,0))
	
	### 6 - Mutations
	par(mar=c(1,0,0,0))
	colAnnotation<-as.matrix(colAnnotation[,toprowOrder])
	if(ncol(colAnnotation)==1){
		colAnnotation<-t(colAnnotation)
	}
	
	plot(0,col="white",xlim=c(0,ncol(colAnnotation)),ylim=c(0,nrow(colAnnotation)),xaxt="n",yaxt="n",type="n",frame.plot=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
	for(i in (nrow(colAnnotation):1)){
		for (j in 1:ncol(colAnnotation)){
			rect(
					xleft=j-0.9,
					ybottom=i-1,
					xright=j+0.1,
					ytop=i,
					col=colAnnotation[nrow(colAnnotation)-i+1,j],
					border="white"
			)
		}
	}
	par(mar=c(0,0,0,0))
	
	
	
	### 7 - Row Dendrogram (or simply row names)
	if(is.null(rowDendro)){
		plot(0,col="white",xlim=c(0,10),ylim=c(0,nrow(plotMatrix)),xaxt="n",yaxt="n",type="n",frame.plot=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
		cex.text<-80/nrow(plotMatrix)
		if((cex.text)>1){cex.text<-1}
		text(rep(10),(1:nrow(plotMatrix)-0.5),labels=rev(rownames(plotMatrix)),pos=2)
	} else {
		plot(as.dendrogram(rowDendro),main="",horiz=TRUE,axes=FALSE,xaxs="i",yaxs="i")
	}
	
	### 8 - Heatmap
	image(t(plotMatrix[nrow(plotMatrix):1,]),col=col,axes=FALSE,breaks=breaks)
}





val2col<-function(z){
	extreme=round(max(abs(z)))
	breaks <- seq(-extreme, extreme, length = nbreaks)
	ncol <- length(breaks) - 1
	col <- colorpanel(ncol,"blue","white","red")
	CUT <- cut(z, breaks=breaks)
	colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
	names(colorlevels)<-names(z)
	return(colorlevels)
}

