kNumStates <- 3
kNumAbnormalStates <- 2

# creates a transition matrix
GetTransitionM = function(distance, homogeneous){
  if (homogeneous){
    prob.normal.normal <- 1 - 1 / 20000
    prob.normal.abnormal <- (1 - prob.normal.normal) / kNumAbnormalStates
    prob.normal.normal <- 1 - prob.normal.abnormal * kNumAbnormalStates
    prob.abnormal.normal <- 1 / 6 
    prob.abnormal.abnormal <- 1 - prob.abnormal.normal
    prob.abnormal.diff.abnormal <- 0
  } else {
    T <- 6
    q <- 1 / T
    p <- 10 ^ -8
    D <- 70000
    f = exp(-distance/D)
    prob.abnormal.abnormal <- f * (1 - q) + (1 - f) * p
    prob.abnormal.normal <- f * q + (1 - f) * (1 - 2 * p)
    prob.abnormal.diff.abnormal <- (1 - f) * p
    prob.normal.normal <- 1 - 2 * p
    prob.normal.abnormal <- p
  }
  transition.probs <- c(prob.abnormal.abnormal, prob.abnormal.normal, 
         prob.abnormal.diff.abnormal, 
         prob.normal.abnormal, prob.normal.normal, prob.normal.abnormal,
         prob.abnormal.diff.abnormal, 
         prob.abnormal.normal, prob.abnormal.abnormal)
 
    transition.m = log(matrix(transition.probs, kNumStates, kNumStates, byrow=TRUE))
    return(transition.m)
}

SumProbabilities <- function(x){
  for (i in 1:(length(x) - 1)) {
    x=sort(x, decreasing = TRUE)
    a = x[1] + log(1+exp(x[2]-x[1]))
    if (length(x) > 2)
      x = c(a, x[(i+2):length(x)]) 
  }     
return(x[1])
}

MultiplyProbabilities <- function(x, y) {x + y}

Forward <- function(emission.probs.mat, distances, homogeneous){
  emission.probs.mat <- as.matrix(emission.probs.mat)
  initial.state.m <- log(rep(0.0075 / kNumAbnormalStates, 3))
  initial.state.m[2] <- log(1 - 0.0075)
  num.windows <- dim(emission.probs.mat)[1]
  forward.m <- matrix(NA, nrow=num.windows, ncol=kNumStates)
  forward.m[1, ] <- initial.state.m + emission.probs.mat[1, ]
  for (i in 2:num.windows){
    temp.m <- forward.m[i-1, ] + GetTransitionM(distances[i], 
                                                homogeneous=FALSE)
    sum.probs <- apply(temp.m, 2, SumProbabilities)
    forward.m[i, ] <- sum.probs + emission.probs.mat[i,] 
  }  
  
  print(SumProbabilities(forward.m[num.windows,]))
  return(forward.m)
}

Backward <- function(emission.probs.mat, distances, homogeneous=FALSE){
  emission.probs.mat <- as.matrix(emission.probs.mat)
  num.windows <- dim(emission.probs.mat)[1]
  initial.state.m <- log(rep(0.0075 / kNumAbnormalStates, 3))
  initial.state.m[2] <- log(1 - 0.0075)
  backward.m <- matrix(NA, nrow = num.windows, ncol=kNumStates)
  
  final.state.m <- rep(0, kNumStates)
  backward.m[num.windows, ] <- final.state.m
  for (i in (num.windows - 1):1){
    temp.m <- GetTransitionM(distances[i+1], FALSE) + 
      matrix(backward.m[i+1, ], 3, 3, byrow=TRUE) +
      matrix(emission.probs.mat[i+1, ], 3, 3, byrow=TRUE)
    backward.m[i, ] <- apply(temp.m, 1, SumProbabilities)  
  }
  final.prob <- backward.m[1, ] + emission.probs.mat[1, ] + initial.state.m
  print(SumProbabilities(final.prob))
  return(backward.m)
}
  
ForwardBackward = function(emission.probs.mat, distances, homogeneous=FALSE){
  
  emission.probs.mat <- as.matrix(emission.probs.mat)
  initial.state.m <- log(rep(0.0075 / kNumAbnormalStates, 3))
  initial.state.m[2] <- log(1 - 0.0075)
  num.windows <- dim(emission.probs.mat)[1]
  forward.m <- matrix(NA, nrow=num.windows, ncol=kNumStates)
  forward.m[1, ] <- initial.state.m + emission.probs.mat[1, ]
  for (i in 2:num.windows){
    temp.m <- forward.m[i-1, ] + GetTransitionM(distances[i], 
                                              homogeneous=FALSE)
    sum.probs <- apply(temp.m, 2, SumProbabilities)
    forward.m[i, ] <- sum.probs + emission.probs.mat[i,] 
  }  
  
  print(SumProbabilities(forward.m[num.windows,]))
  backward.m <- matrix(NA, nrow = num.windows, ncol=kNumStates)
  final.state.m <- rep(0, kNumStates)
  backward.m[num.windows, ] <- final.state.m
  for (i in (num.windows - 1):1){
    temp.m <- GetTransitionM(distances[i+1], FALSE) + 
              matrix(backward.m[i+1, ], 3, 3, byrow=TRUE) +
              matrix(emission.probs.mat[i+1, ], 3, 3, byrow=TRUE)
    backward.m[i, ] <- apply(temp.m, 1, SumProbabilities)  
  }
  browser()
  final.prob <- backward.m[1, ] + emission.probs.mat[1, ] + initial.state.m
  print(SumProbabilities(final.prob))
  
  
} 

GetCNVProbability <- function(forward.m, backward.m, emission.probs.mat, distances, start.target, end.target, state){
  start.forward.prob <- forward.m[start.target, state]
  end.backward.prob <- backward.m[end.target, state]
  segment.prob <- 0
  
  if (end.target > start.target){
    for (i in seq(start.target + 1, end.target)){
      transition.probability <- GetTransitionM(distances[i], homogeneous=FALSE)[state, state]
      emission.probability <- emission.probs.mat[i, state]
      segment.prob <- segment.prob + transition.probability + emission.probability                           
    }
  }
  return(start.forward.prob + end.backward.prob + segment.prob - 
           SumProbabilities(forward.m[dim(emission.probs.mat)[1], ]))
}


Viterbi = function(emission.probs.mat, distances, homogeneous=FALSE){
  emission.probs.mat <- as.matrix(emission.probs.mat)
  num.windows <- dim(emission.probs.mat)[1]
  viterbi.m <- matrix(NA, nrow=num.windows, ncol=kNumStates)
  viterbi.pointers <- matrix(NA, nrow=num.windows, ncol=kNumStates)
  initial.state.m <- log(rep(0.0075 / kNumAbnormalStates, 3))
  initial.state.m[2] <- log(1 - 0.0075)
  dim(initial.state.m) <- c(kNumStates, 1)
  viterbi.m[1, ] <- initial.state.m + emission.probs.mat[1,]
  for (i in seq(2, num.windows)) {
    #browser()
    temp.m <- viterbi.m[i - 1, ] + GetTransitionM(distances[i], homogeneous)
    viterbi.m[i, ] <- apply(temp.m, 2, max)
    emission.probs <- c(emission.probs.mat[i,])
    dim(emission.probs) <- c(kNumStates, 1)
    viterbi.m[i, ] <- viterbi.m[i, ] + emission.probs
    viterbi.pointers[i, ] <- apply(temp.m, 2, which.max)
  }
  viterbi.states = vector(length = num.windows)
  viterbi.states[num.windows] = which.max(viterbi.m[num.windows, ])
  for (i in (num.windows - 1):1) {
    viterbi.states[i] <- viterbi.pointers[i + 1, viterbi.states[i + 1]]
  }
  return(list(viterbi.states=viterbi.states, viterbi.m=viterbi.m))
}


# INPUT
# the first row of counts is:
#chromosome start   end width        gc X1.00079.01 X1.00079.02 X1.00079 X1.00111.01 X1.00111.02
#1          1 69091 70008   918 0.4281046         159         662      636           2          34
#X1.00111 X1.00124.01 X1.00124.02 X1.00124 X1.00148.01 X1.00148.02 X1.00148 X1.00154.01 X1.00154.02
#1       57           1           0      130         218         403      101         300         229
#X1.00154 X1.00159.01 X1.00159.02 X1.00159 X1.00159.1 X1.00170.01 X1.00170.02 X1.00170 X1.00172.01
#1      319         226         413      207        855           0         465      263         224
#X1.00172.02 X1.00222.02 target      
#1         454         282      1 

# test.sample.name="X1.00159.02"

# ref.sample.names =  [1] "X1.00079.01" "X1.00079.02" "X1.00079"    . . . .  (length 24)

callcnvs.Call <- function(counts.normalized.file, fam.file, job.name){
  counts.normalized <- read.table(counts.normalized.file, header=T)
  source("fam.R")
  sample.names <- fam.GetIDs(fam.file)
  source("util.R")
  sample.names.Rformat <- util.GetRFormatFromID(sample.names)
  
  # some of the samples in the fam.file may not be sequenced, so 
  # we select only those that are
  sample.names.Rformat.sequenced <- intersect(sample.names.Rformat, 
                                                 colnames(counts.normalized))
  source("pickrefs.R")
  cluster.data.file <- paste(job.name, ".clusters.txt", sep="")
  pickrefs.WriteClusters(counts.normalized, sample.names.Rformat.sequenced, cluster.data.file)
  clusters <- read.table(cluster.data.file, header=T, 
                         colClasses=c("character", "numeric", "numeric", "integer"))
  clusters$sampleid <- util.GetIDFromRFormat(clusters$sampleid)
  source("fam.R")
  children.ids <- util.GetRFormatFromID(fam.GetChildIDs(fam.file))
  source("xcnv.R")
  xcnv.file <- paste(job.name, ".xcnv", sep="")
  excluded.samples.file <- paste(job.name, ".excluded.samples.txt", sep="")
  xcnv.New(xcnv.file)
  source("emission.R")
  distances <- emission.GetDistances(counts.normalized)
  for (sample.name in head(sample.names.Rformat.sequenced, 9)){
    sample.id <- util.GetIDFromRFormat(sample.name)
    sample.relations <- c(fam.GetMotherID(fam.file, sample.id), fam.GetFatherID(fam.file, sample.id), 
                          fam.GetChildID(fam.file, sample.id))
    sample.cluster <- clusters[clusters$sampleid == sample.id, "cluster"]
    close.samples <- clusters[clusters$cluster == sample.cluster, "sampleid"]
    #print(length(close.samples))
    samples.to.exclude <- c(util.GetRFormatFromID(sample.relations), children.ids, sample.id)
    reference.samples <- setdiff(util.GetRFormatFromID(close.samples), samples.to.exclude)
    if (length(reference.samples) < 8) 
      cat("Excluding sample ", sample.id, "\n", file=excluded.samples.file, append=T)
    else{
      reference.samples.file <- paste(job.name, ".", sample.id, ".reference.samples.txt", sep="")
      var.estimate.file <- paste(job.name, ".", sample.id, ".var.estimate.txt", sep="")
      emission.probs.file <- paste(job.name, ".", sample.id, ".emission.probs.txt", sep="")
      write.table(reference.samples, file=reference.samples.file, quote=F, row.names=F)
      counts.normalized$mean <- apply(counts.normalized[, reference.samples], 1, mean)
      #browser()
      var.estimate <- emission.GetVarEstimate(counts.normalized, reference.samples)
      write.table(var.estimate, file=var.estimate.file, quote=F, row.names=F)
      emission.probs <- GetEmissionProbs(counts.normalized[, sample.name], 
                                       counts.normalized$mean, var.estimate)
      write.table(emission.probs, file=emission.probs.file, quote=T, row.names=F)
      #cat("Emission probabilities calculated\n")
      vlist <- Viterbi(emission.probs, distances)  
      #cat("Viterbi algorithm finished\n")
      viterbi.state <- as.numeric(unlist(vlist["viterbi.states"]))
      #forward.m <- Forward(emission.probs, distances)
      #backward.m <- Backward(emission.probs, distances)
      xcnv.PrintCNVs(sample.name, viterbi.state, counts.normalized, xcnv.file)
    }
  }
  
  
  
#   for (test.name in test.names){
#     cat(as.character(Sys.time()), ": calling CNVs on sample", test.name, "\n")
#     #get emission probabilities
#     emission.probs <- GetEmissionProbs(target.counts, nonzero.ref.counts$mean, var.estimate)
#     
#     distances <- GetDistances(nonzero.ref.counts)
#     #  ForwardBackward(counts.matrix[nonzero.rows, c("del_prob", "normal_prob", "dup_prob")], distances, homogeneous=FALSE)
#     vlist <- Viterbi(emission.probs, distances)  
#     cat("Viterbi algorithm finished\n")
#     viterbi.state <- as.numeric(unlist(vlist["viterbi.states"]))
#     #forward.m <- Forward(emission.probs, distances)
#     #backward.m <- Backward(emission.probs, distances)
#     PrintCNVs(test.name, viterbi.state, nonzero.ref.counts, output.filename)
#   }
#   #counts.matrix[nonzero.rows, c("Vitdelprob", "Vitnormprob", "Vitdupprob")] <- vlist[["viterbi.m"]]
#   #return(GetCNVInfo(counts.matrix[nonzero.rows,]))
#   #browser()  
#   #return(list(nonzero.rows=nonzero.rows, counts.matrix=counts.matrix))
}





PlotTargetData <- function(counts, start.target, end.target, start.cnv, end.cnv, test.sample.name, ref.sample.names){
  par(mfrow=c(1, 1))
  data <- counts[counts$target >= start.target & counts$target <= end.target, ]
  means <- apply(data[, ref.sample.names], 1, mean)
  sd <- sqrt(apply(data[, ref.sample.names], 1, var))
  refs.z.scores <- matrix(NA, dim(data)[1], length(ref.sample.names))
  target.z.score <- numeric(length = dim(data)[1])
  
  for (i in seq(1, dim(data)[1])){
    refs.z.scores[i, ] <- as.numeric((data[i, ref.sample.names] - means[i]) / sd[i])
    target.z.score[i] <- (data[i, test.sample.name] - means[i]) / sd[i]
  }
  ylim <- max(abs(refs.z.scores), abs(target.z.score))
  plot(seq(-6, 6), seq(-6, 6), 
       xlim=c(data[1, "start"], data[dim(data)[1], "start"]), 
       ylim=c(-ylim - 0.1, ylim + 0.1), type="n", xlab="", ylab="Z-score")
  for (i in seq(1, length(ref.sample.names))){
    lines(data[, "start"], refs.z.scores[, i], col="#cdc9c966")
  }
  lines(data[, "start"], target.z.score, col="red")
  points(data[, "start"], rep(-ylim - 0.05, length(data[, "start"])), pch=20)
  cnv.data <- counts[counts$target >= start.cnv & counts$target <= end.cnv, ]
  browser()
  lines( c(cnv.data[1, "start"], cnv.data[dim(cnv.data)[1], "end"]) , c(ylim+0.05, ylim+0.05), lwd=2)
  title(main=paste(counts$chromosome[start.target], ":", data$start[1], "-", data$end[dim(data)[1]], sep=""))
}


