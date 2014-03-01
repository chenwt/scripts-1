##~/tools/R/R_current/bin/Rscript
#-----function
chr2Num = function(chr) {
  if (chr == 'X' || chr =='x') 
  {chr = 23}
  else if (chr == 'Y' || chr == 'y')
  { chr = 24 }
  else if (length(grep(chr,vapply(1:22,as.character,'1'))) > 0 ) 
  { chr = as.integer(chr)}
  else 
  {print("error ")}
  return(chr)
}
blankPlot = function( maxX, maxY){  
  plot(c(-10,0),col="white", 
       xlim = c(0, maxX*1.05), ylim = c(0, maxY+5),
       xaxt="n",yaxt="n",type="n",frame.plot=FALSE,
       xlab="",ylab="", xaxs="i",yaxs="i",mar=c(0, 0, 0, 0))
}
cnv2color = function(str){
  if (str == "DEL") { str = "blue"}
  else if (str == "DUP") {str = "red"}
  else {print("error")}
  return(str)
}


val2col<-function(z, zlim, col = heat.colors(12), breaks){
  #----this function converts a vector of values("z") to a vector of color
  #levels. One must define the number of colors. The limits of the color
  #scale("zlim") or the break points for the color changes("breaks") can 
  #also be defined. when breaks and zlim are defined, breaks overrides zlim.
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  CUT <- cut(z, breaks=breaks)
  colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
  return(colorlevels)
}


##--highlight each rows one by another 
colorRowBg <- function (rowLabels) {
  for (i in 1:length(rowLabels)){
    colBg = ifelse(i %% 2 , "lightblue","white")
    xlim = par()$pin[1]
    rect(0, i-1, xlim, i,border="white",col=colBg)    
  }
}

##--draw lines between column
lineColBorder <- function (coory, i, numPts) {
###-----need further development
  for (i in 1:(length(coory))){ 
    if (i == 1) {
      lines(x=c(coory[i],coory[i]), y = c(0,numPts))
      text(coory[1]/2, numPts + 0.5, labels=names(coory)[1],cex=0.8)
    } else {
      lines(x=c(coory[i],coory[i]), y = c(0,numPts))
      text((coory[i-1] + coory[(i)] )/2, numPts + 0.5, labels=names(coory)[i],cex=0.8)
    }
  }
}
##---axis label
labelTickerXY <- function (rowLabels, colLabels) {
    axis(2,at=0.5:(length(rowLabels)-0.5),cex.axis=0.4, 
         labels = rowLabels,las=2,lwd=0.3, tck = -0.01,  col.ticks="gray")
    axis(3, at=0.5:(length(colLabels)-0.5),cex.axis = 0.7, 
         labels = colLabels,las=2,lwd=0.3, tck = -0.01,  col.ticks="gray")
    axis(1, at= length(colLabels)/2,labels=names(colLabels))    
}

##--plot CNV data
drawHeatmap <- function (data) {
  ##---draw table   
  for ( j in 1:nrow(data)){      
      temp_y = j 
      temp = data[j,]
      for (i in 1:length(temp_y)){
          if ( i == 1) { 
            temp_xs = 0
            temp_xe = i 
          }else{
            temp_xs = i
            temp_xe = i + 1
          }
          temp_col = data[i,j]
          rect(temp_xs , temp_y-1, temp_xe, temp_y, border=temp_col, col=temp_col ) 
      }
  }
}

drawTable <- function (data) {
  ##---draw table   
  for ( j in 1:nrow(data)){      
    temp_y = j 
    temp = data[j,]
    print(temp)
    for (i in 1:length(temp_y)){
      if ( i == 1) { 
        temp_xs = 0
        temp_xe = i 
      }else{
        temp_xs = i
        temp_xe = i + 1
      }
      temp_col = data[i,j]
      print(paste(i, j))
      rect(temp_xs , temp_y-1, temp_xe, temp_y, border="black", col="white" ) 
      text((temp_xs + temp_xe)/2 , temp_y - 0.5, labels=data[i,j])
    }
  }
}


###----plot table
##----load systematic Investor Toolbox
# install.packages("RCurl")
# require(RCurl)
# sit = getURLContent('https://github.com/systematicinvestor/SIT/raw/master/sit.gz', binary=TRUE, followlocation = TRUE, ssl.verifypeer = FALSE)
# con = gzcon(rawConnection(sit, 'rb'))
# source(con)
# close(con)

##---plot data matrix

plot.table = function
(
  plot.matrix,
  smain = '',
  text.cex = 1,
  frame.cell = T,
  highlight = F,
  colorbar = FALSE,
  keep_all.same.cex = FALSE
)
{
  if( is.null(rownames(plot.matrix)) & is.null(colnames(plot.matrix)) ) {
    temp.matrix = plot.matrix
    if( nrow(temp.matrix) == 1 ) temp.matrix = rbind('', temp.matrix)
    if( ncol(temp.matrix) == 1 ) temp.matrix = cbind('', temp.matrix)
    plot.matrix = temp.matrix[-1, -1, drop = FALSE]
    colnames(plot.matrix) = temp.matrix[1, -1]
    rownames(plot.matrix) = temp.matrix[-1, 1]
    smain = temp.matrix[1, 1]
  } else if( is.null(rownames(plot.matrix)) ) {
    temp.matrix = plot.matrix
    if( ncol(plot.matrix) == 1 ) temp.matrix = cbind('', temp.matrix)
    plot.matrix = temp.matrix[, -1, drop = FALSE]
    colnames(plot.matrix) = colnames(temp.matrix)[-1]
    rownames(plot.matrix) = temp.matrix[,1]
    smain = colnames(temp.matrix)[1]
  } else if( is.null(colnames(plot.matrix)) ) {
    temp.matrix = plot.matrix
    if( nrow(temp.matrix) == 1 ) temp.matrix = rbind('', temp.matrix)
    plot.matrix = temp.matrix[-1, , drop = FALSE]
    rownames(plot.matrix) = rownames(temp.matrix)[-1]
    colnames(plot.matrix) = temp.matrix[1, ]
    smain = rownames(temp.matrix)[1]
  }
  plot.matrix[which(trim(plot.matrix) == 'NA')] = ''
  plot.matrix[which(trim(plot.matrix) == 'NA%')] = ''
  plot.matrix[which(is.na(plot.matrix))] = ''
  if(colorbar) {
    plot.matrix = cbind(plot.matrix, '')
    if(!is.null(highlight)) if(!is.logical(highlight)) { highlight = cbind(highlight, NA) }
  }
  nr = nrow(plot.matrix) + 1
  nc = ncol(plot.matrix) + 1
  is_highlight = T
  if(is.logical(highlight)) {
    is_highlight = highlight
    if(highlight) highlight = plot.table.helper.color(plot.matrix)
  }
  if(!is_highlight) {
    plot.matrix.cex = matrix(1, nr = nr, nc = nc )
    plot.matrix_bg.col = matrix('white', nr = nr, nc = nc )
    plot.matrix_bg.col[seq(1, nr, 2), ] = 'lightblue'
    plot.matrix_bg.col[1,] = 'gray';
    plot.table.param( plot.matrix, smain, plot.matrix.cex, plot.matrix_bg.col,
                      frame.cell, keep_all.same.cex)
  } else {
    plot.matrix.cex = matrix(1, nr = nr, nc = nc )
    plot.matrix_bg.col = matrix('white', nr = nr, nc = nc )
    plot.matrix_bg.col[1,] = 'gray'
    plot.matrix_bg.col[2:nr,2:nc] = highlight
    plot.table.param(plot.matrix, smain, plot.matrix.cex, plot.matrix_bg.col,
                     frame.cell, keep_all.same.cex)
  }
  if(colorbar) plot.table.helper.colorbar(plot.matrix);
}
