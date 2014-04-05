def bruteSearch(subseq, seq):
    while subseq[0] in seq:
        index = seq.index(subseq[0])
        if subseq == seq[index:index+len(subseq)]:
            return index
        else:
            seq = seq[index+1:]
    else:
            return -1

def kmpSearch(W, S):
    ''' 
    kmp_search:
    search W against W, return the position of match 
    as the starting of 0
    input:
        an array of characters, S (the text to be searched)
        an array of characters, W (the word sought)
    output:
        an integer (the zero-based position in S at which W is found)
    '''
    define variables:
        an integer, m ← 0 (the beginning of the current match in S)
        an integer, i ← 0 (the position of the current character in W)
        an array of integers, T (the table, computed elsewhere)
    m = 0
    i = 0 
    T = T
    while m + i < len(S):
        if W[i] == S[m + i]:
            if i == len(W) - 1:
                return m
            i = i + 1
        else:
            m = m + i - T[i]
            if T[i] > -1 : 
                i = T[i]
            else :
                i = 0          
    return len(S)

def kmpTable(W):
    '''build kmp_table T for the search algorithm
    input:
        an array of characters, W (the word to be analyzed)
        an array of integers, T (the table to be filled)
    output:
        nothing (but during operation, it populates the table)
    define variables:
        an integer, pos ← 2 (the current position we are computing in T)
        an integer, cnd ← 0 (the zero-based index in W of the next 
character of the current candidate substring)
    (the first few values are fixed but different from what the algorithm 
might suggest)
    '''
    T = [] * len(W); T.append(-1) ; T.append(0)
    while pos < len(W):
        # (first case: the substring continues)
        if W[pos - 1] == W[cnd]: 
            cnd = cnd + 1
            T[pos] = cnd
            pos = pos + 1

        # (second case: it doesn't, but we can fall back)
        elif cnd > 0 :
            let cnd = T[cnd]

        # (third case: we have run out of candidates.  Note cnd = 0)
        else :
            let T[pos] = 0
            pos = pos + 1


