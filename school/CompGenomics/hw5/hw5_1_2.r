##------------------------------
# b.Use (a) to produce several long sequences with k>n, 
# and for compute the k×k matrix that tallies the total frequency counts of 
# pairs of consecutive output symbols along the sequence. Compute eigenvalues for 
# this matrix. What can you empirically say about the eigenvalues? Submit matrices 
# and eigenvalues you have computed

myHW5_2 <- function(l, k, n) {
	inSeq <- myHW5_1(l,k,n)
	# print(inSeq)

	freqMX <- matrix(0,ncol=k,nrow=k)
	rownames(freqMX) <- 1:k
	colnames(freqMX) <- 1:k

	for (i in 2:l) {
		freqMX[inSeq[i-1],inSeq[i]] <- freqMX[inSeq[i-1],inSeq[i]] + 1
	}
	ev <- eigen(freqMX)$values
	
	print(freqMX)
	print(round(ev,0))

}

# l <- 20 # output sequence length
# k <- 5	# observation
# n <- 2  # state
# myHW5_2(l,k,n)
# myHW5_1(l,k,n)