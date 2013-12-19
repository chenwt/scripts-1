getArgs = function(){
	#function used to get commandline argument specified by -- or by sequence
	args = commandArgs(trailingOnly = TRUE)
	hh <- paste(unlist(args),collapse=' ')
	if(length(grep("--",hh)) > 0){
	listOptions <- unlist(strsplit(hh,'--'))[-1]
	optionsArgs <- sapply(listOptions,function(x){
	         unlist(strsplit(x, ' '))[-1]
	        })
	optionsNames <- sapply(listOptions,function(x){
	  option <-  unlist(strsplit(x, ' '))[1]
	})
	names(optionsArgs) <- unlist(optionsNames)
 }else{
 	optionsArgs = args
 }
	return(optionsArgs)
}