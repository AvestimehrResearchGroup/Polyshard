#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Utilities for the variant of the Guruswami-Sudan algorithm using the small
field Kotter-Vardy multiplicity assignment."""

#In this document, using ABC refers to the technical report by Augot, Barbier
#and Couvreur.

from codinglib.listdecoding import *

def gsKV_satisfiable(n,kGRS,tau,q,s):
    """Return whether the given values allows for the construction of a Q
    polynomial. Checks whether there will be more variables than equations.
    """
    gamma = tau/n
    Rr = (kGRS-1)/n
    core = (1-gamma)^2 + gamma^2/(q-1)
    D = s*n*core
    ell = gilt( D/(kGRS-1) )
    return s > 0 and (ell+1)*(D - (kGRS-1)*ell/2) > \
        s*n/2*( (1-gamma)*(s*(1-gamma)+1) + gamma*(s*gamma/(q-1)+1) )


def gsKV_minimal_s(n,kGRS,tau,q):
    """Return a minimal multiplicity for the given parameters"""
    # Keep doubling s to find one that works and then binary
    s = find_minimal_satisfiable(lambda s: gsKV_satisfiable(n,kGRS,tau,q,s))
    assert (gsKV_satisfiable(n,kGRS,tau,q,s))
    return s
        
def gsKV_wu(n,kGRS,tau,q):
    """Use Wu's parameter choice for the multiplicity. This returns s1 and s2,
    where s1 is the multiplicity used for the received field element."""
    d = n-kGRS+1
    tqf = tau^2/(q-1)
    ntS = (n-tau)^2
    tParen = tqf+tau^2+n*d-2*n*tau
    # Equation 34
    C = ligt(
        ( (d*sqrt(tqf+ntS) + sqrt( d^2*(tqf+ntS) - (n^2*(q-1) - n*d*(q-2))*tParen ))
                / 2 / sqrt(2) / tParen
            )^2 - q/8 )
    # Equations 24 and 25
    sqrtCp = sqrt( 2*C + q/4 )
    mDenom = sqrt( (1-tau/n)^2*(q-1)^2 + tau^2/n^2*(q-1))
    m1 = ceil( (1-tau/n)*(q-1)*sqrtCp / mDenom - 1/2)
    m2 = ceil( tau/n*(q-1)*sqrtCp / mDenom - 1/2)
    # Equation 27 and the one two equations above
    ell = floor( ((n-tau)*m1 + tau*m2 -1)/(kGRS-1) )
    #assert (gsKV_satisfiable(n,kGRS,tau,q,s))
    return (m1, m2, ell, C)
                                                 
def gsKV_bernstein(n,d,tau,q):
    """Return the parameters ell, j, k for Bernstein's Alternant decoder"""
    # d is designed minimum distance of the alternant, i.e. t+1 in Bernstein
    # Finds minimal ell which satisfies equation on p. 12 botton
    np = n*(q-1)/q
    Jp = np - sqrt(np*(np-d))
    satisfiable = lambda ell: Jp - tau >= (d/ell + q*n/ell^2)/(2-(tau+Jp)/np)
    ell = find_minimal_satisfiable(satisfiable)
    k = floor(ell*(1-tau/n))
    j = floor(ell*(tau/n)/(q-1))
    return (ell,j,k)

def gsKV_bernstein_satisfiable(n,d,tau,q, ell, j, k):
    """Return whether the parameters for Bernstein's alternant decoder are
    satisfiable"""
    return  (ell >= (q-1)*j + k) \
        and (n*k*(k+1)/2 + n*(q-1)*j*(j+1)/2 + (n-d)*ell*(ell-1)/2
                < ell*(k*(n-tau) + j*tau))

