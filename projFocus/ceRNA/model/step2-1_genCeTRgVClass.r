#/usr/bin/Rscript
#J.HE

#---class constructor
cernet = ""
expmat = ""
samples = 
target   <- setClass("target",
                      slots = c(genename = "character",
                                exp = "numeric",
                                regulators = "character",
                                samples = "character") )

##----construct one ceRNET regulator for one target cancer DEG gene
regulator <- setClass("regulator",
                      slots = c(genename = "character",
                                samples = "character", 
                                exp = "numeric", 
                                cnv = "numeric", 
                                snp = "numeric",
                                som = "numeric",
                                target = "character") )
##---show information of  a regulator: 
setMethod("summary", "regulator",
          function(object){
            cat ("Object of class \"", class(object), "\"\n")
            exp = object@exp
            cnv = object@cnv,
            
          })
##---plot a regulator :
setMethod("plot", "regulator",
          function(object){
            exp = object@exp
          })
  ##exp data

  ##cnv data

  ##snp data

  ##som data

  }
##---initiate a regulator

