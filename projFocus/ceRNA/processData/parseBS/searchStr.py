def bruteSearch(subseq, seq):
    for i in range(len(subseq)):
        if subseq[0] in seq:
            index = seq.index(subseq[0])
            if subseq == seq[index:index+len(subseq)]:
                return index
            else:
                seq = seq[index+1:]
        else:
            return 0 


