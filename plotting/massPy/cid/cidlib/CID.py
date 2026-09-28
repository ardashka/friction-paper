from numpy import log2, mean, random

def cid(hscan, nshuff):
    
    def lz77(sequence):
        sequence = tuple(sequence)
        sub_strings = set()
        ind, inc = 0, 1
        while True:
            if ind + inc > len(sequence):
                break
            sub_str = sequence[ind : ind + inc]
            if sub_str in sub_strings:
                inc += 1
            else:
                sub_strings.add(sub_str)
                ind += inc
                inc = 1
        return len(sub_strings)
    
    def cid(hscan):
        C, L = lz77(hscan), len(hscan)
        return C*(log2(C) + 2*log2(L/C)) / L
    
    def cid_shuff(hscan):
        shuffles = []
        rng = random.default_rng()
        for _ in range(nshuff):
            rng.shuffle(hscan)
            shuffles.append( cid(hscan) )
        return mean(shuffles)
    
    ccid = cid(hscan)
    
    return ccid, 1 - ccid / cid_shuff(hscan, nshuff)


def cid_corr(hscan, nshuff):
    """ WARNING :: This mehtod hasn't been tested yet! """
    qs = []
    while len(hscan) > 1:
        qs.append( cid(hscan, nshuff) )
        del hscan[::-2]
    return qs
