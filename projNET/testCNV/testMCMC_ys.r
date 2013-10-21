#!/bin/Rscirpts
#Author: Jing He
#Date: 
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


a1 = args[1]
b1 = args[2]

a2 = arge[3]
b2 = arge[4] 

n = 50
paiA = c(1:20) * 0.05

est = c(1:length(paiA)) 

llnull = est
llalter = llnull

i = 1
for (pai in paiA) 
{
	a = rbeta(n * pai, a1, b1)
	b = rbeta(n*(1 -pai), a2, b2)
	
	d = c(a,b)
	
	pt = c(0:100)*0.01
	ha = rep(0,length(pt)) # log likelihood
	j = 1
	
	llnull[i] = 0
	for (x in d) {
		llnull[i] = llnull[i] + log(dbeta(x, a2, b2))
	}

	
	llalter[i] = 1
	for (p in pt) {
		h = 0 
		for (x in d) {
			lh = p * dbeta(x, a1, b1) + (1-p)*dbeta(x, a2-1, b2-1)
			h = h + log(lh)
			llalter[i]  = llalter[i] + lh
		}
			
		ha[j] = h 		
		j = j + 1
	}
	llalter[i] = log( llalter[i] / length(pt))
	temp = ha[1]
	est[i] = pt[1]
	for (j in 2:length(pt)) {
		if (ha[j] > temp) {
			temp = ha[j]
			est[i] =pt[j] 
		}
		
	}

	i = i + 1
}