#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Implementation of decoding algorithms Reed-Solomon up to half the
minimum distance"""

from codinglib.util import *
from codinglib.code import *
from codinglib.module import *
from codinglib.bma import bma

def forney_algorithm(C, r, Lambda, Omega):
    """Correct the received word `r` given the error locator `Lambda` and the
    error evaluator `Omega` using the Forney algorithm."""
    roots = [ i for i in range(C.n) if Lambda(1/C.alphas[i]).is_zero() ]
    if not len(roots) == Lambda.degree():
        raise DecodingFailedError("Lambda is not an error locator")
    dLambda = Lambda.derivative()
    vs = C.parity_colmults()
    c = [ ri for ri in r ] # copy r
    for i in roots:
        e = -(C.alphas[i] / vs[i]) * (Omega(1/C.alphas[i]) / dLambda(1/C.alphas[i]))
        c[i] = c[i] - e
    # OBS: Has not been checked whether it is a codeword yet.
    return vector(c)

class GRSDecoderGao(DecoderBlockCodeAbstract):
    r"""
    Implementation of the Gao decoding algorithm for GRS codes, from the paper
      Gao, Shuhong. "A New Algorithm for Decoding Reed-Solomon Codes."
      In Communications, Information and Network Security, 55-68. Springer, 2003.
    As opposed to that paper, however, the implementation uses the
    language of modules and weak Popov form (following jsrn phd
    thesis), instead of the language of Extended Euclidean
    algorithm.

    INPUT:

    - ``C`` - The code to decode
    - ``reencode=False`` - Optional. Specifies whether to apply reencoding before solving the key equation.
    """

    def __init__(self, C, reencode=False):
        self.C = C
        PF.<x> = C.F[]
        self.PF = PF
        self.reencode = reencode
        if reencode:
            # Precompute stuff to make reencoding faster on the fly
            self.G_r = prod(x-alpha for alpha in C.alphas[C.k:])
            self.g = prod(x-alpha for alpha in C.alphas[:C.k])
            hs_r = [ self.G_r//(x-a) for a in C.alphas[C.k:] ] 
            # Lagrangian base functions
            self.Hs_r = [ hs_r[i-C.k]/(hs_r[i-C.k](C.alphas[i]))/self.g(C.alphas[i]) for i in range(C.k,C.n) ]
        else:
            self.G = prod(x-alpha for alpha in C.alphas)

    def decode_to_information(self, r):
        if self.reencode:
            return self._decode_reencode(r)
        else:
            return self._decode_raw(r)

    def _decode_raw(self, r):
        C = self.C
        R = self.PF.lagrange_polynomial([(C.alphas[i], r[i]/C.colmults[i]) for i in range(C.n)])
        M = matrix(self.PF, 2, 2, [ 1, R, 0, self.G])
        ws = [C.k, 0]
        module_apply_weights(M, ws)
        module_weak_popov(M)
        i = 0 if LP(M.row(0))==0 else 1
        module_remove_weights(M, [C.k,0])
        L,Lf = M.row(i)
        if L.divides(Lf):
            f = self.PF(Lf/L)
            return vector(f.list())
        else:
            raise DecodingFailedError

    def _decode_reencode(self, r):
        C = self.C
        f_r = self.PF.lagrange_polynomial([ (C.alphas[i], r[i]/C.colmults[i]) for i in range(0, C.k) ])
        r_r = vector([ C.F.zero() ]*C.k + [ r[i]/C.colmults[i] - f_r(C.alphas[i]) for i in range(C.k,C.n) ])

        R_r = sum( r_r[i] * self.Hs_r[i-C.k] for i in range(C.k,C.n) )
        x = self.PF.gen()
        Mw = matrix(self.PF, 2, 2, [[1, R_r],[0, self.G_r]])
        module_weak_popov(Mw)
        j = 0 if LP(Mw.row(0))==0 else 1
        (L,T) = Mw.row(j)  #If L=Lambda then T = Lambda*f_f/g
        Tt = T*self.g
        if L.divides(Tt):
            f_f = Tt//L
            return poly2list(f_f + f_r, C.k)
        else:
            raise DecodingFailedError

    #TODO: Make static
    def is_list_decoder(self):
        return False

    def decoding_radius(self):
        return (self.C.d-1)//2

    def __str__(self):
        return "Gao-decoder up to d/2 for %s" % self.C




class GRSDecoderSyndrome(DecoderBlockCodeAbstract):
    r"""
    Implementation of the Syndrome key equation decoding algorithm
    for GRS codes, originally described by Berlekamp:
       Berlekamp, Elwyn R. Algebraic Coding Theory. Aegean Park Press, 1968.

    INPUT:

    - ``C`` - The code to decode
    - ``solver="bma"`` - Optional, either "bma" or "ea". Decides which
    algorithm to use for solving the key equation.
    """

    def __init__(self, C, solver='bma'):
        self.C = C
        # TODO: Check that C is a (subfield subcode of) a GRS code
        assert not 0 in C.alphas \
            , "Syndrome decoding is not implemented for codes with zero as an evaluation point"
        PF.<x> = C.F[]
        self.PF = PF
        self.mod = x^(C.n-C.k)
        self.solver = solver
        C.parity_colmults() # cache-calculate parity column multipliers

    def decode(self, r):
        if self.solver=="bma":
            return self._decode_bma(r)
        else:
            return self._decode_ea(r)

    def _decode_ea(self, r):
        C = self.C
        S = C.syndrome(r, direct=False)
        M = matrix(self.PF, 2, 2, [ 1, S, 0, self.mod])
        module_weak_popov(M)
        i = 0 if LP(M.row(0))==0 else 1
        Lambda,Omega = M.row(i)
        return forney_algorithm(self.C, r, Lambda, Omega)

    def _decode_bma(self, r):
        C = self.C
        S = C.syndrome(r, direct=False)
        (Lambda_seq, L, _) = bma(S)
        Lambda = self.PF(Lambda_seq)
        Omega = (self.PF(S) * Lambda) % self.mod
        return forney_algorithm(self.C, r, Lambda, Omega)

    def is_list_decoder(self):
        return False

    def decoding_radius(self):
        return (self.C.d-1)//2

    def __str__(self):
        return "Syndrome key equation decoder up to d/2 for %s" % self.C
