#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

r"""The Extended Euclidean algorithm with stopping conditions.
Note that using modules over F[x] and weak Popov form is a more flexible
approach."""

import codinglib.util

def ea(a, b, t=0, dt=None): 
    """Euclid's Algorithm for polynomials. t gives optionally the degree of the
    remainder for which to break, and dt gives optionally the difference in
    degrees of remainder and the intermediate v poly for which to continue."""
    # The polynomials each iteration satisfies
    # sc = a*uc + b*vc
    (sc, sl) = (b, a)
    (uc, ul) = (0, 1)
    (vc, vl) = (sc.parent()(1), 0)
    i = 1
    if dt is None:
        cont = lambda: sc.degree() >= t
    else:
        cont = lambda: sc.degree() - vc.degree() >= dt
    while cont():
        (q, sn) = sl.quo_rem(sc)
        (sc,sl) = (sn, sc)
        (uc, ul) = (ul - q*uc , uc)
        (vc, vl) = (vl - q*vc , vc)
        i=i+1
        #print "s[%d] = %s\nq[%d] = %s\nu[%d] = %s\nv[%d] = %s" \
        #             % (i,sc.degree(),  i,q.degree(),  i,uc.degree(),
        #                     i,vc.degree())
    return ((sc,uc,vc), (sl, ul, vl))

def ea_equiv(d, syndromes): 
    """Euclids Algorithm for RS-decoding (Sugiyama). Mostly implemented as
    described in Equivalence paper.
    Input:  The code's minimum distance and the syndromes
    Output: Error locator, Error evaluator and Error co-evaluator"""
    PF.<x> = syndromes[0].parent()[] # Init F[x] with x as formal param
    t = d//2
    A = Matrix(2, 2, 1) # 2x2 identity
    p_old = x^(2*t)
    p = list2poly(syndromes,x)
    while p.degree() >= t:
        q = p_old//p
        (p_old, p) = (p, p_old-q*p)
        A = Matrix( ((0, 1), (1, -q)) )*A
    scale = A[1,1](0)
    Lambda = A[1,1]/scale
    Omega  = p/scale
    Psi    = A[1,0]/scale
    return (Lambda, Omega, Psi)
