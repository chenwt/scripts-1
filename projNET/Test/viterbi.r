
apply(ndata[1:5,3:4],1,function(x){
				choose(x[1],x[2]) * exp((lbeta(x[2]+alpha,x[1]-x[2]+beta)-lbeta(alpha,beta)))
				#lbeta to tackle underflow
				})


lapply(states,FUN=function(x){start[x] * tranMX[,x] * emiMX[1,x]})