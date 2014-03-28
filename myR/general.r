
##---get command line arguments

getArgs = function(){
  #base function used to get commandline argument specified by -- or by sequence
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


jxy = function(...){
  ss = unlist(list(...))
  temp = ""
  res = paste(ss, collapse="")
  return(res)
}


setRootd = function(){
  ###used to setup rootDir when working on mapped disk/remote server
  sysInfo = Sys.info()
  if(sysInfo['sysname']=="Darwin" ){
    print("working from MacOS")
    rootd = "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  }else if(sysInfo['sysnameß']=="Linux" ){
    print("working from Linux")
    rootd = "/ifs/home/c2b2/ac_lab/jh3283/projFocus/"
  }
  return(rootd)
}
