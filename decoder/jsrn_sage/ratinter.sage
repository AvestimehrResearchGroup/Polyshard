#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Rational interpolation using the Wu algorithm"""
from codinglib.util import *
from codinglib.listdecoding import *

#####################################################
### CHOOSING PARAMETERS
#####################################################


def rat_satisfiable(n,w,tau,s,l):
    """Returns true iff the given s and l satisfies the constraints for the
    given n,w and tau"""
    return s > 0 and l > 0 and s*(s+1)*n/2 < (l+1)*(s*tau - l*w/2)

def rat_params(n,w,tau):
    """Return my closed expressions for satisfiable s and l, from the Goppa
    paper"""
    #if not list_decodable(n,n-k+1,tau):
    #    raise ValueError(
    #        "tau does not satisfy requirements wrt. n and k:\n%i < tau < %i" %
    #            list_decoding_range(n,d))
    if w==0:
        s = 1
        l = ligt(n/tau-1)
    else:
        smin = w*(n-tau)/(tau^2-n*w)
        s = ligt(smin)
        eps = s - smin
        D = eps*(tau^2 - n*w)*s + w^2/4
        l = floor( tau/w*s + 1/2 - sqrt(D)/w )
    #print (n,k,tau,s,l,w)
    assert rat_satisfiable(n,w,tau,s,l)
    return (s,l)

def rat_params_old(n, w, tau):
    """Calculate the s and l parameters given n, w and tau, using the closed
    expressions of Wu's original paper (and my thesis)"""
    #if not list_decodable(n,n-k+1,tau):
    #    raise ValueError(
    #        "tau does not satisfy requirements wrt. n and k:\n%i < tau < %i" %
    #            list_decoding_range(n,d))
    if w==0:
        s = 1
        l = ligt(n/tau-1)
    else:
        s = floor(tau*(tau-w)/(tau^2 - n*w))
        l = ceil(s*tau/w - 1)
    #print (n,k,tau,s,l,w)
    assert rat_satisfiable(n,w,tau,s,l)
    return (s,l)

def rat_params_min_list(n,w,tau):
    """Find the lowest list size and a corresponding multiplicity for the given
    parameters."""
    # It is probably near the l given by the equations, so start here
    (firsts,firstl) = rat_params(n,w,tau)
    def try_l(l):
        (mins, maxs) = solve2deg_int(n , n-2*tau*(l+1) , l*(l+1)*w)
        if mins <= l and maxs > 0:
            return max(1, mins) 
        else:
            return None
    l = find_minimal_satisfiable(try_l, firstl)
    s = try_l(l)
    if not l is None and not s is None:
        #print (n,k,tau,s,l,w)
        assert rat_satisfiable(n,w,tau,s,l)
        return (s, l)
    else:
        return None

def rat_params_min_mult(n,w,tau):
    """Find the lowest multiplicity and a corresponding list size for the given
    parameters."""
    # It is probably near the s given by the equations, so start here
    (firsts,firstl) = rat_params(n,w,tau)
    def try_s(s):
        (minl, maxl) = solve2deg_int( w,  w-2*s*tau,  n*s*(s+1)-2*s*tau )
        if maxl >= minl and minl >= s:
            return max(minl, s) 
        else:
            return None
    s = find_minimal_satisfiable(try_s, firsts)
    assert rat_satisfiable(n,w,tau, s,try_s(s))
    return (s, try_s(s))

def rat_params_trifonov(n,w,tau):
    """Return Trivonov's expressions for the multiplicity and list size. Note
    that these are sometimes erroneous, e.g. (n,k,tau) = (64,4,46), as the
    interval specified for l does not contain an integer"""
    #if not list_decodable(n,n-k+1,tau):
    #    raise ValueError(
    #        "tau does not satisfy requirements wrt. n and k:\n%i < tau < %i" %
    #            list_decoding_range(n,d))
    if w==0:
        return (1, ligt(n/tau-1))
    s = ligt( w*(n-tau + sqrt( n^2 - 2*tau*n + w*n )) / (2*(tau^2-w*n)) )
    D = (w+2*s*tau)^2 - 4*w*n*s*(s+1.)
    l = ligt( (2*s*tau - w - sqrt(D)) /2/w )
    return (s,l)

def rat_system_size(n, w, tau, (s,l)):
    """Return the size of the Guruswami-Sudan system for the given
    parameters: no. of equations, no. of variables"""
    noeqs = n*s*(s+1)/2
    novars = s*tau*(l+1) - l*(l+1)*w/2
    return (noeqs, novars)

def rat_stats(n, w, tau):
    """Collect some statistics on the rational interpolation problem served.
    Returns a dict with the following entries:
        mins, minl : minimal mult and list size
        exprs, exprl : mult and list size when using the closed expressions
        noeqs: no. of equations to solve to construct Q (using exprs, exprl)
        novars: no. of monomials when constructing Q, using the formula
        maxtrixsize: no. of cells in the matrix to solve for finding Q
    """
    stats = dict()
    (mins, minl) = rat_params_min_mult(n,w,tau)
    stats['mins'] = mins
    stats['minl'] = minl
    (exprs, exprl)= rat_params(n,w,tau)
    stats['exprs'] = exprs
    stats['exprl'] = exprl
    (stats['noeqs'], stats['novars']) = rat_system_size(n,w,tau,(exprs,exprl))
    stats['matrixsize'] = stats['noeqs']*stats['novars']
    return stats


#####################################################
### CONSTRUCTING AN INTERPOLATION POLYNOMIAL
#####################################################

def rat_monomial_list(maxdeg, l, (wy, wz)):
    """Return a list of the (x,y,z) powers of all partially homogenized
    monomials in F[x,y,z] whose (1,wy,wz)-weighted degree is less than maxdeg
    and whose (y+z)-degree = l."""
    mons = []
    for y in range(0, l+1):
        z = l-y
        for x in range(0,  ceil(maxdeg - y*wy - z*wz)):
            mons.append((x, y, z))
    return mons

def rat_interpol_matrix(points, s, mons):
    """Return the interpolation matrix whose nullspace gives the coefficients
    for all interpolation polynomials with the given parameters. The ith column
    will be the coefficients on the ith monomial in mons."""
    n = len(points)
    def eqs_affine(x0, theta0, flip=False):
        """Make equation for the affine point x0, theta0. Return a list of
        dictionary from monomial (as tuple of powers) to coefficient, where the
        monomials are all in mons. Each dictionary corresponds to one equation."""
        eqs = []
        for i in range(0, s):
            for j in range(0, s-i):
                eq = dict()
                for mon in mons:
                    ihat = mon[0]
                    jhat = mon[1] if not flip else mon[2]
                    if ihat >= i and jhat >= j:
                        icoeff = binomial(ihat, i)*x0^(ihat-i) \
                                    if ihat > i else 1
                        jcoeff = binomial(jhat, j)*(theta0^(jhat-j)) \
                                    if jhat > j else 1
                        eq[mon] = jcoeff*icoeff
                eqs.append(eq)
        #print "Eqs ", eqs
        return eqs
    def eqs_projective((x,y,z)):
        """Return a list of equations for the projective point. Each equation
        is representated a list of coefficients whose dot-product with mons
        must equal zero."""
        if not z.is_zero():
            aeqs = eqs_affine(x, y/z)
        else:
            aeqs = eqs_affine(x, 0, flip=True)
        return [ [aeq.get(mon, 0) for mon in mons] for aeq in aeqs ]
    return flatten_once([ eqs_projective(point) for point in points ])

def rat_construct_Q(points, tau, (s,l), weights):
    """Calculate the partially homogeneous interpolation polynomial Q(x,y,z)
    for the parameters given. points is a list of triples (xi,yi,zi) such that
    Q(xi,yi,zi) = 0 with multiplicity s"""
    F = points[0][0].parent()
    mons = rat_monomial_list(tau*s, l, weights)
    desnMons = s*tau*(l+1) - l*(l+1)*(weights[0]+weights[1])/2
    #print "tau*s = %d, (s, l)=(%s,%s),  (wy,wz)=%s, mons=%s (%s), eqs=%s" \
    #       % (s*tau, s, l, weights, len(mons), desnMons, len(points)*s*(s+1)/2)
    #Below: Faster but less readable coeff finder. Should convert well to
    #       Cython.
    #THERE IS A BUG IN THIS SOMEWHERE!!!!
    #nmons = len(mons)
    #n = len(points)
    #eqs_point = int(s*(s+1)/2)
    #neqs = n*eqs_point
    #print "Matrix dimensions: %s x %s" % (neqs, nmons)
    #splitieq = list(flatten_once(
    #    [ [(x,y) for x in range(0,s-y)] for y in range(0,s) ]))
    #def get_monom_coeff(ipoint, ieq, mon):
    #    (x0,y0,z0) = points[ipoint]
    #    i, j = splitieq[ieq]
    #    ihat = mon[0]
    #    if z0.is_zero():
    #        jhat = mon[2]
    #        theta = F.zero()
    #    else:
    #        jhat = mon[1]
    #        theta = y0/z0
    #    if ihat >= i and jhat >= i:
    #        icoeff = binomial(ihat, i)*x0^(ihat-i) if ihat > i else 1
    #        jcoeff = binomial(jhat, j)*(theta^(jhat-j)) if jhat > j else 1
    #        return jcoeff*icoeff
    #    return 0
    #M = matrix(neqs, nmons, lambda r,c: \
    #        get_monom_coeff(int(r)/eqs_point, int(r) % eqs_point, mons[c]))
    M = matrix(list(rat_interpol_matrix(points, s, mons)))
    #print "Constructed matrix (%s x %s)" % (M.nrows(), M.ncols())
    Sp = M.right_kernel()
    sol = Sp.random_element()
    while sol.is_zero():
        print "rat_construct_Q: Found zero as random element. Trying again."
        # Picking out e.g. element 1 directly seems to run into an infinite
        # loop for large matrices.
        sol = Sp.random_element()
    # Construct the Q polynomial
    PF.<x,y,z> = F[]
    Q = sum([ x^mons[i][0]*y^mons[i][1]*z^mons[i][2]*sol[i]
                for i in range(0, len(mons)) ])
    return Q

def rat_find_factors(Q,w1,w2):
    """Returns a list of tuples (f1, f2) where (y*f1 - z*f2) | Q and deg(f1) <=
    w1 and deg(f2) <= w2"""
    Px.<x> = Q.base_ring()[]
    Pyz.<y,z> = Px[]
    facs = Q.factor(proof=False) # fingers crossed
    goodFacs = []
    for [fac,mult] in facs:
        yzfac = Pyz(fac)
        if set([y,z]).issuperset(yzfac.monomials()):
            funy = Px(yzfac.coefficient([1,0])) # coeffs are in F[x]
            funz = Px(yzfac.coefficient([0,1]))
            if funy.degree() <= w1 and funz.degree() <= w2:
                goodFacs.append((funy, funz))
    return goodFacs
