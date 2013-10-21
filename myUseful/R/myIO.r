myWTablet <- function(inobj,outfile){
	write.table(inobj,"outfile",quote=FALSE,col.name=FALSE,row.name=FALSE,sep="\t" )
}

myWTablec <- function(inobj,outfile){
	write.table(inobj,"outfile",quote=FALSE,col.name=FALSE,row.name=FALSE,sep="," )
}

# myRTable <- read.table()
