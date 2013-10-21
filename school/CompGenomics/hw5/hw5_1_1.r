##------------------------------
# a.	Write a program that receives input integers l,k,n, 
# chooses random n×n transition and n×k emission matrices 
# and produces an HMM output sequence 
# of length l of symbols 1,…,k .



myHW5_1 <- function(l, k, n) {
	# ##------------------------------parameter
	# l <- 10 # output sequence length
	# k <- 2	# observation
	# n <- 2  # state
	##------------------------------transition matrix
	set.seed(100)
	tranMX <- matrix(runif(n^2),ncol=n)
	tranMX <- tranMX / rowSums(tranMX)
	colnames(tranMX) <- 1:n
	rownames(tranMX) <- 1:n 
	##------------------------------emission matrix
	set.seed(100)
	emisMX <- matrix(runif(n*k),ncol=k)
	emisMX <- emisMX / rowSums(emisMX)
	colnames(emisMX) <- 1:k
	rownames(emisMX) <- 1:n
	##------------------------------start
	start <- runif(n)
	start <- start/sum(start)

	outMX  <- matrix(NA, nrow=l,ncol=n)

	outState  <- vector(length=l)
	outObserve <- vector(length=l)

	##------------------------------initiate
	temp <- log(start) + log(emisMX)
	# rownames(outMX) <- 1:l
	# colnames(outMX) <- 1:n
	# outState[1] <- which.max(temp[,which.max(colSums(temp))])
	outState[1] <- sample(1:n,1)
	# outObserve[1] <- which.max(colSums(temp))
	outObserve[1] <- sample(1:k,1)
	outMX[1,] <- temp[,outState[1]]

	for(i in 2:l) {
		# print("previous out")
		# print(outMX[i-1,])
		temp1 <-  outMX[i-1,] + log(tranMX) 
		temp <- colSums(temp1) + log(emisMX)

		# outMX[i,] <- temp + log(apply(emisMX,1,max))
		# outState[i] <- which.max(temp[,which.max(colSums(temp))])
		outState[i] <- sample(seq(1,n),1)
		outMX[i,] <- temp[,outState[i]]
		# outObserve[i] <- which.max(colSums(temp))
		outObserve[i] <-sample(seq(1,k),1)
	}
	return(outObserve)
}

