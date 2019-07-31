#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

r"""The Berlekamp-Massey algorithm."""

def bma(seq):
    """Berlekamp-Massey's Algorithm, Blahut 2008, p. 157, slightly modified
    Input:  A sequence seq over field (length = 2t)
            var: the string name of the variable to use in the polynomials
    Output: (A, B): A is the shortest sequence such that the convolution of A
            over seq at all shifts is zero. B is the companion polynomial
            computed by the BMA. Note that this is the reverse of the
            characteristic polynomial's coefficients

    EXAMPLES
        sage: bma([1,2,1,2,1,2])
        [1,0,-1]
        sage: bma([GF(7)(1),19,1,19])
        [1,0,6]
        sage: bma([2,2,1,2,1,191,393,132])
        [1, -36727/11711, 34213/5019, 7024942/35133, -335813/1673]
    """
    F = seq[0].parent()
    zero, one = F.zero(), F.one()
    A = [one]+[zero]*len(seq)
    B = [one]+[zero]*len(seq)
    t = [zero]*(len(seq)+1)
    L = 0
    shift = 1
    for r in range(0, len(seq)):
        disc = zero
        for i in range(0,L+1):
            disc += A[i]*seq[r-i]
        if disc.is_zero(): 
            shift += 1
        else:
            if L >= (r+1)-L:
                for i in range(0, (r+1)-L+1-shift):
                    A[i+shift] -= disc*B[i]
                shift += 1
            else:
                for i in range(0, shift):
                    t[i] = A[i]
                    A[i] /= disc
                for i in range(0, (r+1)-L+1-shift):
                    t[i+shift] = A[i+shift] - disc*B[i]
                    A[i+shift] /= disc
                L = r+1-L
                A,B,t = t,A,B
                shift=1

    #Post processing so users don't get confused
    A = A[:L+1]
    B = [zero]*(shift-1) + B[:len(seq)-L+1-(shift-1)]
    return A, L, B

def bma_poly(S):
    """BMA working on polynomials: S is a polynomial and returns A of least
    degree such that deg(A*S % x^{len(S)}) < deg(A). B is the companion
    polynomial calculated by the BMA."""
    Aseq, L, Bseq = bma(S.list())
    PF = S.parent()
    return PF(Aseq), PF(Bseq)

def produces(seq, A):
    """True iff the linear recurrence produces the input sequence"""
    if isinstance(A, tuple):
        # Guess that its (polynomial, length)
        (Ap, L) = A
        A = Ap.list() + [0]*(L-Ap.degree())
    L = len(A)
    A = A[:] #copy
    A.reverse()
    good = True
    for i in range(L, len(seq)):
        if not sum([ a*b for (a,b) in zip(A, seq[i-L:i+1])]).is_zero():
            print "Failed %i'th cipher (%i)" % (i+1, seq[i])
            return False
    return True
