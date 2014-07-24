rm(list=ls())
require(cgdsr)

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

# Get list of cancer studies at server
mycgid = "brca"
mycdIndx = grep(mycgid, getCancerStudies(mycgds)[,1])


# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[11,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[1,1]

# Get data slices for a specified list of genes, genetic profile and case list
fcon =  file("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/062014/wgsVar/tcag_pub_gene.list",open='r')
genelist = readLines(fcon)
close(fcon)
genelist = unlist(strsplit(genelist,"\n"))

mut = na.omit(getProfileData(mycgds,genelist[1:100], mygeneticprofile,mycaselist)

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)

# documentation
help('cgdsr')
help('CGDS')



# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

# Get list of cancer studies at server
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[2,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[4,1]

# Get data slices for a specified list of genes, genetic profile and case list
getProfileData(mycgds,c('BRCA1','BRCA2'),mygeneticprofile,mycaselist)

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)