import timeit
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

# def searchUTR(subseq, seq):
#     tlen = len(seq)
#     slen = len(subseq)
#     i = 0 
#     while i < tlen - slen:
#         if not (seq[i].islower() and seq[i+slen].islower() ):
#             i = i + 1
#             continue
#         else:
#             if subseq == seq[i: i+slen]:
#                 yield i
#             i = i + 1

def searchUTR(subseq, seq):
    subseq = subseq.lower()
    tlen = len(seq)
    slen = len(subseq)
    i = 0 
    while i < tlen - slen:
        seq[i:i+slen].lower() 
        if subseq == seq[i: i+slen]:
            yield i
        i = i + 1


