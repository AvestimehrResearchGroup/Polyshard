#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Functionality for Power decoding of Reed-Solomon codes.
Power decoding was called "Decoding by virtual extension to an interleaved
code" in the initial description of [SSB].

In the initial algorithm of [SSB], the decoding radius is bounded by (roughly)
the Sudan radius. In [N15], it was shown how to incorporate multiplicities to
decode up to the Johnson radius. The implementation here supports that.

REFERENCES:
 - Schmidt, G., V.R. Sidorenko, and M. Bossert. "Syndrome Decoding of
   Reed-Solomon Codes Beyond Half the Minimum Distance Based on Shift-Register
   Synthesis." IEEE Transactions on Information Theory 56, no. 10 (2010):
   5245-52. doi:10.1109/TIT.2010.2060130.

 - Nielsen, Johan S. R. "Power Decoding Reed-Solomon Up to the Johnson
   Radius.". Submitted to IEEE Transactions of Information Theory, 2015.
"""

from codinglib.util import *
from codinglib.code import DecoderBlockCodeAbstract, DecodingFailedError
from codinglib.key2d import KeyEquation2D
from codinglib.module import *
from codinglib.rs_mindist import forney_algorithm


### Not-exported functions

def power_decoding_radius_from_parameters_real(n, k, (s,ell)):
    return (2*ell-s+1)/2/(ell+1)*n - ell/2/s*(k-1) - ell/(ell+1)

def power_decoding_radius_from_parameters(n,k, (s,ell)):
    return floor(power_decoding_radius_from_parameters_real(n, k, (s,ell)))

def make_powers(p, max_pow):
    ps = [ p.parent().one(), p ]
    for i in range(1,max_pow):
        ps.append(ps[i]*p)
    return ps

### Decoder classes

class GRSDecoderPower(DecoderBlockCodeAbstract):
    """
    Power Decoder for GRS codes, Gao style.
    Supports multiplicities.
    """

    def __init__(self, C, params):
        self.C = C
        self.s = params[0]
        self.ell = params[1]
        Px.<x> = C.F[]
        self._Px = Px

    def precompute(self):
        x = self._Px.gen()
        G = prod(x-alpha for alpha in self.C.alphas)
        self._Gs = make_powers(G, self.s)

    @staticmethod
    def decoding_radius_from_parameters_real(C, (s,ell)):
        return power_decoding_radius_from_parameters_real(C.n, C.k, (s,ell))

    @staticmethod
    def decoding_radius_from_parameters(C, (s,ell)):
        return power_decoding_radius_from_parameters(C.n, C.k, (s,ell))

    def decoding_radius_real (self):
        return GRSDecoderPower.decoding_radius_from_parameters_real(self.C, (self.s, self.ell))

    def decoding_radius(self):
        return GRSDecoderPower.decoding_radius_from_parameters(self.C, (self.s, self.ell))

    def _build_module(self, r):
        """Build the punctured module which is a basis of the solution module"""
        C = self.C
        s,ell = self.s, self.ell
        tau = self.decoding_radius()

        Gs = self._Gs
        R = self._Px.lagrange_polynomial(zip(C.alphas, [ r[i]/C.colmults[i] for i in range(C.n) ]))
        Rs = make_powers(R, ell)
        def cells(i, j):
            if i==0 and j==0:
                return 1
            if i < s:
                t = j
                if i <= t:
                    return (binomial(t, i)*Rs[t-i]*Gs[i]) % Gs[s]
                else:
                    return 0
            if i == j:
                return Gs[s]
        M = matrix(self._Px, ell+1, ell+1, cells)
        #print poly_degs(M)
        # print
        return M

    def _short_vector(self, M):
        """Find the shortest weighted vector in the module"""
        k = self.C.k; ell= self.ell
        ws = [ell*(k-1) + 1] + [ (ell-t)*(k-1) for t in range(1,ell+1) ]
        # print ws
        module_weak_popov(M, weights=ws)
        # print poly_degs(M)
        # print 
        sol_row = module_row_with_LP(M, 0, weights=ws)
        return M.row(sol_row)


    def decode_to_information(self, r):
        if not hasattr(self, '_Gs'):
            self.precompute()
        M = self._build_module(r)
        sol =  self._short_vector(M)
        lambdaS = sol[0]
        phi = sol[1]
        # print "polys: ", lambdaS, "\n", phi
        if lambdaS.is_zero():
            raise Exception("Internal error: Found zero error locator")
        if lambdaS.divides(phi):
            return poly2list(phi//lambdaS, self.C.k)
        else:
            raise DecodingFailedError("Power Decoding Failed")

    def is_key_equation_solution(self, r, lambdas, psis, debug=0):
        """Return whether the given lambdas and psis solve the key equations
        corresponding to the received word r"""
        C = self.C
        R = self._Px.lagrange_polynomial([ (C.alphas[i], r[i]/C.colmults[i]) for i in range(C.n) ])
        Rs = make_powers(R, s)
        s,ell = self.s, self.ell
        for t in range(1, ell+1):
            lsum = sum( lambdas[i] * binomial(t,i) * Rs[t-i] * self._Gs[i] for i in range(min(s-1, t)+1) )
            if t < s:
                if not lsum == psis[t-1]:
                    if debug:
                        print "Key equation %s (equality) was not satisfied" % t
                        if (self._Gs[s]).divides(lsum - psis[t-1]):
                            print "\tBut congruence WAS satisfied"
                    return False
            else:
                if not (self._Gs[s]).divides(lsum - psis[t-1]):
                    if debug:
                        print "Key equation %s (equality) was not satisfied" % t
                    return False
        return True



    def is_list_decoder(self):
        return False




class GRSDecoderPowerSyndrome(DecoderBlockCodeAbstract):
    """
    Power Decoder for GRS codes, Syndrome style.
    Supports multiplicities.
    """

    def __init__(self, C, params):
        self.C = C
        if C.F.zero() in C.alphas:
            raise ValueError("Syndrome decoding only works for GRS codes not using 0 as an evaluation point")
        self.s = params[0]
        self.ell = params[1]
        Px.<x> = C.F[]
        self._Px = Px

    def precompute(self):
        x = self._Px.gen()
        s = self.s
        G = prod(x-alpha for alpha in self.C.alphas)
        self._Gs = make_powers(G, s)
        rG = G.reverse()
        tau = self.decoding_radius()
        rGi = rG.inverse_mod(x^(s*(self.C.n-self.C.k)+1))
        rGis = make_powers(rGi, s)
        #rGis is the inverse of (the reverse of G) modulo x^(large enough power)
        self._rGis = rGis

    @staticmethod
    def decoding_radius_from_parameters_real(C, (s,ell)):
        return power_decoding_radius_from_parameters_real(C.n, C.k, (s,ell))

    @staticmethod
    def decoding_radius_from_parameters(C, (s,ell)):
        return power_decoding_radius_from_parameters(C.n, C.k, (s,ell))

    def decoding_radius_real (self):
        return GRSDecoderPower.decoding_radius_from_parameters_real(self.C, (self.s, self.ell))

    def decoding_radius(self):
        return GRSDecoderPower.decoding_radius_from_parameters(self.C, (self.s, self.ell))

    def _build_module(self, r):
        """Build the punctured module which is a basis of the solution module"""
        C = self.C
        n,k = C.n, C.k
        s,ell = self.s, self.ell
        tau = self.decoding_radius()

        x = self._Px.gen()
        Gs = self._Gs
        rGis = self._rGis
        R = self._Px.lagrange_polynomial(zip(C.alphas, [ r[i]/C.colmults[i] for i in range(n) ]))
        Rs = make_powers(R, ell)

        mods = [0]+[ x^(t*(n-k)) for t in range(1,s+1) ] \
                  +[ x^(s*n - t*(k-1) - 1) for t in range(s+1, ell+1) ]
        # print mods

        def cells(i, j):
            if i==j and i<s:
                return 1
            if i < s and j < s:
                return 0
            if i < s and j >= s:
                t = j-s+1
                if i <= t:
                    T = (Rs[t-i] % Gs[s-i]) * Gs[i]
                    if t <= s:
                        rev_deg = (t-i)*(n-1) + i*n
                        xpow = 0
                    else:
                        rev_deg = s*n-1
                        xpow = i
                    # print "i=%s, t=%s,  T=%s" % (i, t, T.degree())
                    rT = T.reverse(degree=rev_deg)
                    return (binomial(t, i)*x^xpow*rT*rGis[min(s,t)]) % mods[t] 
                else:
                    return 0
            if i == j:
                # print i, len(mods)
                return mods[i-s+1]
        M = matrix(self._Px, ell+s, ell+s, cells)
        # print poly_degs(M)
        # print
        return M

    def _short_vector(self, M):
        """Find the shortest weighted vector in the module"""
        n, k = self.C.n, self.C.k;
        s, ell= self.s, self.ell
        tau = self.decoding_radius()
        target_degs = [s*tau - i for i in range(s)] + [0]*(s-1) + [ s*tau - s ] + [ s*tau - 1 for t in range(s+1, ell+1) ]
        max_deg = max(target_degs)
        ws = [ max_deg - tdeg for tdeg in target_degs ]
        ws[0] += 1
        # print ws
        # print
        module_weak_popov(M, weights=ws)
        # print poly_degs(M)
        sol_row = module_row_with_LP(M, 0, weights=ws)
        # print "Row ", sol_row
        return M.row(sol_row)

    def decode(self, r):
        if not hasattr(self, '_Gs'):
            self.precompute()
        Px = self._Px
        M = self._build_module(r)
        sol =  self._short_vector(M)
        lambdaS = sol[0].reverse()
        phi = sol[1].reverse()
        if lambdaS.is_zero() or phi.is_zero():
            raise DecodingFailedError("Found polynomials are not non-zero")
        div = lambdaS.gcd(phi)
        # print div
        Lambda, Omega = Px(lambdaS/div) , Px(phi/div)

        if not Lambda.degree() == lambdaS.degree()/self.s:
            raise DecodingFailedError("Found polynomial is not the s'th power of an error locator")

        roots = lambdaS.roots(multiplicities=False)
        if len(roots) != Lambda.degree():
            raise DecodingFailedError("Found polynomial is not an error locator")

        c = forney_algorithm(self.C, r, Lambda, Omega)
        if self.C.iscodeword(c):
            raise DecodingFailedError("The Forney algorithm did not correct into a codeword")
        return c

        #TODO
        # f = self._Px.lagrange_polynomial([ (C.alphas[i], 

        # print poly_degs(M)
        # print sol
        # print "polys: ", lambdaS, "\n", phi
        # if lambdaS.is_zero():
        #     raise Exception("Internal error: Found zero error locator")
        # if lambdaS.divides(phi):
        #     return poly2list(phi//lambdaS, self.C.k)
        # else:
        #     raise DecodingFailedError("Power Decoding Failed")

    # def is_key_equation_solution(self, r, lambdas, psis, debug=0):
    #     """Return whether the given lambdas and psis solve the key equations
    #     corresponding to the received word r"""
    #     C = self.C
    #     R = self._Px.lagrange_polynomial([ (C.alphas[i], r[i]/C.colmults[i]) for i in range(C.n) ])
    #     s,ell = self.s, self.ell
    #     for t in range(1, ell+1):
    #         lsum = sum( lambdas[i] * binomial(t,i) * R^(t-i) * self._Gs[i] for i in range(min(s-1, t)+1) )
    #         if t < s:
    #             if not lsum == psis[t-1]:
    #                 if debug:
    #                     print "Key equation %s (equality) was not satisfied" % t
    #                     if (self._Gs[s]).divides(lsum - psis[t-1]):
    #                         print "\tBut congruence WAS satisfied"
    #                 return False
    #         else:
    #             if not (self._Gs[s]).divides(lsum - psis[t-1]):
    #                 if debug:
    #                     print "Key equation %s (equality) was not satisfied" % t
    #                 return False
    #     return True



    def is_list_decoder(self):
        return False


    



#####
#Some decoding radius/parameter functions (instead of actually constructing a real decoder)
#####
#TODO: These are somewhat deprecated and do not support multiplicities.

def power_correction_radius(n,k,ell=None):
    """Return the error correction radius for Power decoding (slightly less than
    Sudan) for the parameters. If ell is not given, return the highest possible
    correction radius"""
    if not ell:
        ell = power_best_list_size(n,k)
    return floor( ell/(ell+1)*n - ell/2*(k-1) - ell/(ell+1) )

def power_best_list_size(n,k):
    """Return the list size which results in the greatest error correction
    radius for Power decoding"""
    #This is the equation from Sidorenko and Schmidt.
    #return floor( (sqrt((k+3)^2 + 8*(k-1)*(n-1)) - (k+3))/(2*(k-1)) )
    #This is the equation from jsrn phd, p. 97
    r = 1/2 + 1/(k-1)
    return floor( sqrt(r^2 + 2*(n-2)/(k-1)) - r)

def power_virtual_received(r, C, t):
    """Return the t'th power of the normalised received word with GRS base code C."""
    vs = C.colmults
    return vector( power_discr(r[i]/vs[i], t) for i in range(0,C.n) )


#
### Gao variant of Power decoding
#

def power_gao_KE(C, r, ell):
    """Return the Key Equation corresponding to Power Gao decoding"""
    PF.<x> = C.F[]
    G = prod(x-alpha for alpha in C.alphas)
    Rs = [ PF.lagrange_polynomial(zip(C.alphas, power_virtual_received(r,C,t))) for t in range(1, ell+1) ]
    return KeyEquation2D(sigma=ell, rho=1, Gs=[G]*ell, Ss=[Rs], nu=1, etas=[ell*(C.k-1)+1], ws=[ (ell-j)*(C.k-1) for j in range(1,ell+1) ])


#
### Syndrome variant of Power decoding
#

def power_virtual_syndrome_len(C, t):
    """Return the length of the t'th virtual syndrome. Note, this is the number
    of coefficients in that virtual syndrome polynomial and NOT its degree
    (which is one less)"""
    return C.n-t*(C.k-1)-1

def power_virtual_syndrome(r, C, t):
    """Return the t'th virtual syndrome with base code C"""
    slen = power_virtual_syndrome_len(C, t)
    n = C.n
    alphas = C.alphas
    etas = C.multipliers_product()
    rt = power_virtual_received(r, C, t)
    PF = PolynomialRing(r[0].base_ring(), 'x')
    return PF([ sum( rt[i]*etas[i]*alphas[i]^j for i in range(0, n) ) for j in range(0, slen) ])
    
def power_syndrome_KE(C, r, ell):
    """Return the Key Equation corresponding to Power Syndrome decoding"""
    PF.<x> = C.F[]
    mods = [ x^power_virtual_syndrome_len(C, t) for t in range(1, ell+1) ]
    Ss = [ power_virtual_syndrome(r, C, t) for t in range(1, ell+1) ]
    return codinglib.key2d.KeyEquation2D(sigma=ell, rho=1, Gs=mods, Ss=[Ss], nu=1, etas=[0], ws=[0]*ell)
