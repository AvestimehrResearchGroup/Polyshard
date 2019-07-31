#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Utility functions related to list decoding"""
from codinglib.util import *

def list_decoding_range(n, d, q=None):
    """Return the minimal and maximal number of errors correctable by a
    Johnson-distance list decoder beyond half the minimal distance."""
    if q is None:
        return ( ligt((d-1)/2) , gilt(n - sqrt(n*(n-d))) )
    else:
        f = (q-1.)/q
        return ( ligt((d-1)/2) , gilt(f*(n-sqrt(n*(n-d/f)))) )

def list_decodable(n,d,tau):
    """Return whether tau is an acceptable list decoding radius for the given n
    and minimum distance"""
    (mint,maxt) = list_decoding_range(n,d)
    return tau >= mint and tau <= maxt

def received_word_with_list(C, tau, minWeight=None, minTau=None, minList=2, zeroIsInList=True):
    """Find a received word which results in lists of size at least minList
    when decoded to radius tau in the code C. Returns None if no such word
    exists. If minTau is given, all words will be at least minTau from the
    received word. If minWeight is given and minList=2, the non-zero word found
    with low weight will have at least minWeight.
    The algorithm always initially finds a word which is close to the zero
    word. If zeroIsInList is False then a random codeword from C will be added
    to the received word.
    NOTE: Exponential in C's dimension as the entire code is tabulated!"""
    from random import shuffle
    if minWeight and minTau and minWeight < 2*minTau:
        minWeight = 2*minTau
    if (minList==2):
        zero = C.F.zero()
        if hasattr(C, 'low_weight_codeword'):
            # Use C's own low weight codeword function to find the words
            if minWeight:
                lowWord = C.low_weight_codeword(minWeight)
            elif minTau:
                lowWord = C.low_weight_codeword(2*minTau)
            else:
                lowWord = C.low_weight_codeword()
        else:
            # Find a low weight codeword brute force by tabulating C
            if minTau:
                check = lambda(v): weight(v) >= 2*minTau and weight(v) <= 2*tau
            else:
                check = lambda(v): nonzero_word(v) and weight(v) <= 2*tau
            clst = C.list()
            shuffle(clst)
            lowWord = (v for v in clst if check(v)).next()
        lowWordWeight = weight(lowWord)
        if not minTau:
            minTau = 0
        midWordWeight = randint( max(minTau, lowWordWeight-tau),
                                 min(tau, lowWordWeight-minTau))
        midWord = subset_word(lowWord, midWordWeight)
        c1, c2 = vector(C.F, len(lowWord)) , lowWord
        if not zeroIsInList:
            c = C.random_codeword()
            midWord = midWord + c
            c1, c2 = c1+c, c2+c
        return (midWord, [c1, c2])
    else:

        raise NotImplementedError("Lists of size >2 not automatically found")
        #smallWords = filter(lambda v: weight(v) <= 2*tau, allWords);
        #if len(smallWords) >= minList:
