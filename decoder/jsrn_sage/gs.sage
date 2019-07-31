#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""The Guruswami-Sudan algorithm for Reed-Solomon codes: implementation and utilities"""
from codinglib.util import *
from codinglib.code import *
import codinglib.listdecoding
import codinglib.rootfinding
import codinglib.module
import codinglib.codeTesting


#TODO: Restructure construction functions to return a Q in F[x][y]
#TODO: Rename l to ell everywhere

MES_IMPOSSIBLE_PARAMS = "Impossible parameters for the Guruswami-Sudan algorithm"

def gs_satisfiable(n,k,tau,s,l):
    """Returns whether the given parameters satisfy the governing equation of
    Guruswami-Sudan (e.g. jsrn msc thesis p. 36)"""
    return l > 0 and s > 0 and n*s*(s+1) < (l+1)*(2*s*(n-tau) - (k-1)*l)

def gs_params_old(n,k,tau):
    """"Calculate values of (s,l) like those given in the paper of
    Guruswami-Sudan"""
    if k==1:
        return (1,1)
    kk = k-1
    t = n-tau
    s = 1 + floor( (kk*n + sqrt(kk^2*n^2 + 4*(t^2-kk*n))) / (2*(t^2 - kk*n)) )
    l = floor((s*t - 1)/kk)
    assert gs_satisfiable(n,k,tau,s,l) , MES_IMPOSSIBLE_PARAMS
    return (s, l)

def gs_params(n,k,tau):
    """Calculate values of (s,l) like those given in jsrn phd thesis."""
    w = k-1
    atau = n-tau
    smin = tau * w / (atau^2 - n*w)
    s = floor(1+smin)
    D = (s-smin)*(atau^2-n*w)*s + w^2/4 
    l = floor( atau / w * s + 1/2 - sqrt(D)/w )
    assert gs_satisfiable(n,k,tau,s,l) , MES_IMPOSSIBLE_PARAMS
    return (s, l)

def gs_params_mceliece(n,k,tau):
    """Calculate the values of (s,l) like those given by McEliece in 2003
    technical report, using the expression for s isolated by Cassuto and
    Bruck.
    I'm pretty sure there's a mistake in McEliece's multiplicity expression. If
    one plugs the found multiplicity and list size into gs_satisfiable, it
    often fails. McEliece uses a pretty different setup so it is hard to see
    where his requirements come from, but it would be very surprising if it
    could be better than the usual GSA derivation. An example which fails, and
    where (44) of McEliece is true but probably shouldn't be is
    (n,k,tau)=(1024,512,260) where McEliece finds mult 2."""
    kk = k-1
    tt = tau+1
    #s = ligt( kk*(tt + sqrt(n*(2*tt+kk-n)))/2/((n-tt)^2-kk*n) )
    s = ligt( kk*(tau + sqrt(n*(2*tau+kk-n)))/2/((n-tau)^2-kk*n) )
    #(s,dummy) = gs_minimal_mult(n,k,tau)
    l = floor(  sqrt(n/kk*s*(s+1)+((kk+2)/2/kk)^2) - (kk+2)/2/kk )
    #assert gs_satisfiable(n,k,tau,s,l) # This often fails
    return (s, l)

def gs_params_list_leeosullivan(n, k, s):
    """Lee and O'Sullivan's list size from 'List Decoding of Reed-Solomon Codes
from a Groebner Basis Perspective.' 2008 """
    N = n * s * (s+1)/2 + 1
    return gilt( sqrt(2*N/(k-1) + 1/4) - 1/2 )

def gs_params_leeosullivan(n, k, tau):
    """Lee and O'Sullivan's list size with our multiplicity."""
    (s,our_ell) = gs_params(n,k,tau)
    l = gs_params_list_leeosullivan(n,k,s)
    assert gs_satisfiable(n,k,tau,s,l) , MES_IMPOSSIBLE_PARAMS
    return (s, l)

def gs_minimal_list(n,k,tau):
    """Calculate the minimal list size for the parameters along with a
    satisfiable multiplicity. Returns (s, l). This is slow (linear in the list
    size)."""
    # It is probably near the l given by the equations, so start here
    (firsts,firstl) = gs_params(n,k,tau)
    def try_l(l):
        (mins,maxs) = solve2deg_int(n, n-2*(l+1)*(n-tau), (k-1)*l*(l+1))
        if maxs > 0 and maxs >= mins:
            return max(1, mins) 
        else:
            return None
    l = find_minimal_satisfiable(try_l, firstl)
    s = try_l(l)
    assert gs_satisfiable(n,k,tau,s,l) , MES_IMPOSSIBLE_PARAMS
    return (s, l)

def gs_minimal_mult(n,k,tau):
    """Calculate the minimal mult size for the parameters along with a
    satisfiable list size. Returns (l, s). This is slow (linear in the list
    size)."""
    # It is probably near the s given by the equations, so start here
    (firsts,firstl) = gs_params(n,k,tau)
    def try_s(s):
        (minl,maxl) = \
                solve2deg_int(k-1 , k-1-2*s*(n-tau) , n*s*(s+1)-2*s*(n-tau) )
        if maxl >= 0 and maxl >= minl:
            return max(1, minl) 
        else:
            return None
    s = find_minimal_satisfiable(try_s, firsts)
    l = try_s(s)
    assert gs_satisfiable(n,k,tau,s,l) , MES_IMPOSSIBLE_PARAMS
    return (s,l)

def gs_decoding_radius(n, k, l=None, s=None):
    """Return the maximal decoding radius of the Guruswami-Sudan decoder and
    minimal parameters attaining this. If ell is set, the maximal radius using
    this ell; same for s."""
    def get_tau(s,l):
        # jsrn msc thesis, eqn. 3.3.1, page 36
        if s<=0 or l<=0:
            return -1
        return gilt( n - n/2*(s+1)/(l+1) - (k-1)/2*l/s )
    if l==None and s==None:
        tau = codinglib.listdecoding.list_decoding_range(n,n-k+1)[1]
        return (tau, gs_minimal_list(n,k,tau))
    if l!=None and s!=None:
        return (get_tau(s,l), (s,l))
    if s!= None:
        # maximising tau under condition
        # n*(s+1 choose 2) < (ell+1)*s*(n-tau) - (ell+1 choose 2)*(k-1)
        # knowing n and s, we can just minimise
        # ( n*(s+1 choose 2) + (ell+1 choose 2)*(k-1) )/(ell+1)
        # Differentiating and setting to zero yields ell best choice:
        lmax = sqrt(n*s*(s+1)/(k-1)) - 1
        #the best integral value will be
        (l,tau) = find_integral_max(lmax, lambda l: get_tau(s,l))
        assert gs_satisfiable(n,k,tau,s,l) , MES_IMPOSSIBLE_PARAMS
        #Note that we have not proven that this ell is minimial in integral
        #sense! It just seems that this most often happens
        return (tau,(s,l))
    if l!= None:
        # Acquired similarly to when restricting s
        smax = sqrt((k-1)/n*l*(l+1))
        (s,tau) = find_integral_max(smax, lambda s: get_tau(s,l))
        assert gs_satisfiable(n,k,tau,s,l) , MES_IMPOSSIBLE_PARAMS
        return (get_tau(s,l), (s,l))

def gs_system_size(n, k, tau, (s,l)):
    """Return the size of the Guruswami-Sudan system for the given
    parameters: no. of equations, no. of variables"""
    noeqs = n*s*(s+1)/2
    novars = (l+1)*(s*(n-tau) - (k-1)*l/2)
    return (noeqs, novars)

def gs_stats(n, k, tau):
    """Collect some statistics on the polynomial interpolation problem served.
    Returns a dict with the following entries:
        mins, minl : minimal mult and list size
        exprs, exprl : mult and list size when using the closed expressions
        noeqs: no. of equations to solve to construct Q using exprs, exprl
        nomineqs: no. of equations to solve to construct Q using mins, minl
        novars: no. of monomials when constructing Q, using the formula
        maxtrixsize: no. of cells in the matrix to solve for finding Q with
                     exprs, exprl
        minmaxtrixsize: no. of cells in the matrix to solve for finding Q with
                     mins, minl
    """
    stats = dict()
    (mins, minl) = gs_minimal_list(n,k,tau)
    stats['mins'] = mins
    stats['minl'] = minl
    (exprs, exprl)= gs_params(n,k,tau)
    stats['exprs'] = exprs
    stats['exprl'] = exprl
    (stats['noeqs'], stats['novars']) = gs_system_size(n,k,tau,(exprs,exprl))
    (stats['minnoeqs'], stats['minnovars']) = gs_system_size(n,k,tau,(mins,minl))
    stats['matrixsize'] = stats['noeqs']*stats['novars']
    stats['minmatrixsize'] = stats['minnoeqs']*stats['minnovars']
    return stats


#####################################################
### CONSTRUCTING AN INTERPOLATION POLYNOMIAL
#####################################################
#
# The Q construction functions all have the same signature such that
# they can replace each other seamlessly in any setting, e.g. in
# gs_grs_decode, also in this file:
#   points * tau * (s * l) * weight_y  --> Q-poly
#
# TODO: The use of tau here is only used by the linalg constructor for
# computing how many linear restrictions shoul be used. It is a
# negligible micro-optimisation for the case where one decodes fewer
# errors than possible: one can instead compute the least possible Q
# from the input parameters, implying the maximal tau. This is
# implicitly done by the other constructors.
#
# TODO: Currently some of the Q-constructors output Q in F[x,y] while
# others output in F[x][y]. They should all use the latter.
#

### Construction by solving a linear system of equations
def gs_monomial_list(maxdeg, l, wy):
    """Return a list of the (x,y) powers of all monomials in F[x,y] whose
    (1,wy)-weighted degree is less than maxdeg and whose y-degree <= l."""
    mons = []
    for y in range(0, l+1):
        for x in range(0,  ceil(maxdeg - y*wy)):
            mons.append((x, y))
    return mons

def gs_interpol_matrix_by_mons(points, s, mons):
    """Return the interpolation matrix whose nullspace gives the coefficients
    for all interpolation polynomials, given the list of monomials allowed. The
    ith column will be the coefficients on the ith monomial in mons."""
    n = len(points)
    def eqs_affine(x0,y0):
        """Make equation for the affine point x0, y0. Return a list of
        equations, each equation being a list of coefficients corresponding to
        the monomials in mons."""
        eqs = []
        for i in range(0, s):
            for j in range(0, s-i):
                eq = dict()
                for mon in mons:
                    ihat = mon[0]
                    jhat = mon[1]
                    if ihat >= i and jhat >= j:
                        icoeff = binomial(ihat, i)*x0^(ihat-i) \
                                    if ihat > i else 1
                        jcoeff = binomial(jhat, j)*(y0^(jhat-j)) \
                                    if jhat > j else 1
                        eq[mon] = jcoeff*icoeff
                eqs.append([eq.get(mon, 0) for mon in mons])
        return eqs
    return flatten_once([ eqs_affine(*point) for point in points ])

def gs_interpol_matrix_problem(points, tau, (s,l), wy):
    """Return the linear system of equations which Q should be a solution to.
    Returns a matrix M and a list of monomials mons, where a vector in the
    right nullspace of M corresponds to an interpolation polynomial $Q$, by the
    $i$'th element being the coefficient of the $i$'th monomial in mons of
    $Q$."""
    mons = gs_monomial_list((len(points)-tau)*s, l, wy)
    M = matrix(list(gs_interpol_matrix_by_mons(points, s, mons)))
    return (M, mons)

def gs_construct_Q_from_matrix(M, mons):
    """Given the interpolation matrix problem and the corresponding list of
    monomials, return a satisfactory Q polynomial."""
    if M.nrows() >= M.ncols():
        raise Exception("More rows than columns! Bailing")
    Sp = M.right_kernel()
    sol = Sp.an_element()
    #TODO: Option to pick out minimal element?
    while sol.is_zero():
        print "rat_construct_Q: Found zero as random element. Trying again."
        # Picking out e.g. element 1 directly seems to run into an infinite
        # loop for large matrices.
        sol = Sp.random_element()
    # Construct the Q polynomial
    PF.<x,y> = M.base_ring()[]
    Q = sum([ x^mons[i][0]*y^mons[i][1]*sol[i] for i in range(0, len(mons)) ])
    return Q

def gs_construct_Q_linalg(points, tau, (s,l), wy):
    """Calculate an interpolation polynomial Q(x,y) for the parameters
    given by solving a linear system of equations.
    points is a list of tuples (xi,yi) such that Q(xi,yi) = 0 with multiplicity
    s.
    Shorthand for calling
    gs_construct_Q_from_matrix(*gs_interpol_matrix_problem(...))
    """
    return gs_construct_Q_from_matrix(
                *gs_interpol_matrix_problem(points, tau, (s,l), wy))



    
### Construction by module minimising Lee and O'Sullivan's basis
def gs_lee_osullivan_module(points, tau, (s,l), wy):
    """Return the analytically straight-forward basis for the module containing
    all interpolation polynomials, as according to Lee and O'Sullivan"""
    F = points[0][0].parent()
    PF.<x> = F[]
    R = PF.lagrange_polynomial(points)
    G = prod( x - points[i][0] for i in range(0, len(points)) )
    PFy.<y> = PF[]
    ybasis = [ (y-R)^i*G^(s-i) for i in range(0, s+1) ] \
            + [ y^(i-s)*(y-R)^s for i in range(s+1, l+1) ]
    def pad(lst):
        return lst + [0]*(l+1-len(lst))
    modbasis = [ pad(yb.coefficients(sparse=False)) for yb in ybasis ]
    return Matrix(PF, modbasis)

def gs_construct_Q_lee_osullivan(points, tau, (s,l), wy, minimiser=codinglib.module.module_weak_popov):
    """Module minimise the Lee-O'Sullivan module in order to find a satisfactory Q."""
    F = points[0][0].parent()
    M = gs_lee_osullivan_module(points, tau, (s,l), wy)
    weights = [ i*wy for i in range(0,l+1) ]
    codinglib.module.module_apply_weights(M, weights)
    minimiser(M)
    # Construct Q as the element of the row with the lowest weighted degree
    degs = [ (i,LT(M.row(i)).degree()) for i in range(0,l+1) ]
    best = min(degs, key=lambda (i,d): d)[0]
    codinglib.module.module_remove_weights(M, weights)
    Qlist = M.row(best)
    PFxy.<xx,yy> = F['x,y']
    Q = sum( yy^i*PFxy(Qlist[i]) for i in range(0,l+1) )
    return Q



### Construction by Koetter--Nielsen--Hoeholdt
def gs_construct_Q_knh(points, tau, (s,l), wy):
    """Construct interpolation polynomial Q by building a Groebner
    basis of all interpolation polynomials bottom-up using the
    Koetter--Nielsen--Hoeholdt algorithm.
    This implementation uses the Hasse matrix speed-up described in
      Fast Koetter--Nielsen--Hoeholdt Interpolation in the Guruswami--Sudan Algorithm
      Johan S. R. Nielsen
      ACCT 2014
    Note: This is NOT the fast algorithm of that paper, but simply the
    regular KNH with the Hasse matrix speed-up.
    """
    F = points[0][0].parent()
    PF.<x> = F[]
    PFy.<y> = PF[]
    basis = [ y^i for i in range(l+1) ]
    derivatives = [ (i,j) for i in range(0,s) for j in range(0,s-i) ] # order is important
    for point in points:
        hasse = [ gs_Q_hasse_derivatives(basis[i], point, s) for i in range(l+1) ]
        for deriv in derivatives:
            evals = [ hasse[i][deriv] for i in range(l+1) ]
            #TODO: Cache and update weighted degrees?
            degs = [ gs_Q_weighted_degree(basis[i], wy) for i in range(l+1) ]
            t = -1
            for i in range(l+1):
                if not evals[i].is_zero() and (t==-1 or degs[i] < degs[t]):
                    t = i
            if t > -1:
                # fix all other indices with basis[t]
                for i in range(l+1):
                    if i != t and not evals[i].is_zero():
                        coeff =  -evals[i]/evals[t]
                        basis[i] += basis[t]*coeff
                        hasse[i] += hasse[t]*coeff
                # fix basis[t] and hasse[t]
                basis[t] *= x-point[0]
                hasse[t] = zero_matrix(F, 1, s).stack(hasse[t].delete_rows([s-1]))
    degs = [ (i, gs_Q_weighted_degree(basis[i],wy)) for i in range(l+1) ]
    Q = basis[min(degs, key=lambda (i, deg): deg)[0]]
    return Q


def gs_construct_Q_knh_orig(points, tau, (s,l), wy):
    """Construct interpolation polynomial Q by building a Groebner
    basis of all interpolation polynomials bottom-up using the
    Koetter--Nielsen--Hoeholdt algorithm.
    This implementation DOES NOT use the Hasse matrix speed-up, and so
    is quite a lot slower than gs_construct_Q_knh (seemingly on all
    parameters)
    """
    F = points[0][0].parent()
    # We represent Q-polys as list of lists (2d arrays), where there
    # are always l+1 outer lists, each of them representing an F[x]
    # polynomial.
    def normalize(px):
        """Cut off leading 0s in a polynomial in list representation"""
        while px and px[-1].is_zero():
            px.pop()
    def add_polys(p1, p2, coeff):
        """Compute p1 + coeff*p2 on 2d array representations"""
        p = []
        for yi in range(len(p1)):
            py1 = p1[yi]; py2 = p2[yi]
            py = []
            for xi in range(min(len(py1),len(py2))):
                py.append(py1[xi] + coeff*py2[xi])
            for xi in range(len(py2), len(py1)):
                py.append(py1[xi])
            for xi in range(len(py1), len(py2)):
                py.append(coeff*py2[xi])
            if len(py1) == len(py2):
                normalize(py)
            p.append(py)
        return p
    def shift_poly(p, x0):
        """Compute p*(x-x0) on 2d array representation"""
        q = []
        for yi in range(len(p)):
            py=p[yi]
            qy=[]
            if py:
                qy.append(-x0*py[0])
                for xi in range(1,len(py)+1):
                    c = py[xi-1] - x0*py[xi] if xi<len(py) else py[xi-1]
                    qy.append(c)
            q.append(qy)
        return q
    def weighted_degree(p):
        deg = -1
        for yi in range(len(p)):
            t = yi*wy + len(p[yi]) if p[yi] else -1
            if t > deg:
                deg = t
        return deg
    basis = [ [ [] for j in range(i) ] + [[F.one()]] + [[] for j in range(i+1,l+1)] for i in range(l+1) ]
    derivatives = [ (i,j) for i in range(0,s) for j in range(0,s-i) ] # order is important
    # Cache binomial coefficients since this is crazy slow in Sage
    maxx = s*len(points)
    binomials = [ [ binomial(a,b) for b in range(s) ] for a in range(maxx) ]
    for (x0,y0) in points:
        #Cache powers of x0 and y0, mainly to avoid function calls
        xpows = [ F.one(), x0 ]+[ x0^i for i in range(2, maxx) ]
        for (dy,dx) in derivatives:
            # Compute Hasse derivatives
            ypows = [ F.one(), y0 ]+[ y0^i for i in range(2, l+1) ]
            evals = []
            for b in basis:
                d = F.zero()
                for yi in range(dy, len(b)):
                    dt = F.zero()
                    for xi in range(dx, len(b[yi])):
                        dt += binomials[xi][dx] * xpows[xi-dx] * b[yi][xi]
                    dt *= binomials[yi][dy] * ypows[yi-dy]
                    d += dt
                evals.append(d)
            degs = [ weighted_degree(basis[i]) for i in range(l+1) ]
            t = -1
            for i in range(l+1):
                if not evals[i].is_zero() and (t==-1 or degs[i] < degs[t]):
                    t = i
            if t > -1:
                # fix all other indices with basis[t]
                for i in range(l+1):
                    if i != t and not evals[i].is_zero():
                        coeff =  -evals[i]/evals[t]
                        basis[i] = add_polys(basis[i], basis[t], coeff)
                # fix basis[t]
                basis[t] = shift_poly(basis[t],x0)
    degs = [ (i, weighted_degree(basis[i])) for i in range(l+1) ]
    Q = basis[min(degs, key=lambda (i, deg): deg)[0]]
    PF.<x> = F[]
    PFy.<y> = PF[]
    return PFy([ PF(qi) for qi in Q ])





### Construction by Fast variant of Koetter--Nielsen--Hoeholdt
def gs_knh_fast_mod_key(points):
    """Returns the mod-tree lookup key for a given point list, for the precomputed KNH mod tree"""
    # assumes all points have different x coordinate
    return (points[0][0], len(points))

def gs_knh_fast_mod_tree(points, s, x):
    """Pre-compute the tree of moduli used for the KNH-fast algorithm.
    Returns a dictionary mapping (x-coord * length) to moduli, where
    the key indicates the first x-coordinate and its length of a given
    point-list encountered in the tree."""
    lib = dict()
    def fill_points(points): #mirror structure of gs_knh_fast_tree
        if len(points)==1:
            lib[gs_knh_fast_mod_key(points)] = (x-points[0][0])^s
        else:
            half = len(points)//2
            first, second = points[:half] , points[half:]
            fill_points(first)
            fill_points(second)
            lib[gs_knh_fast_mod_key(points)] = lib[gs_knh_fast_mod_key(first)] * lib[gs_knh_fast_mod_key(second)]
    fill_points(points)
    return lib
        
def gs_knh_fast_reduce_basis(basis, points, aux):
    modulo =  aux['mods'][gs_knh_fast_mod_key(points)]
    return basis.apply_map(lambda q: q % modulo, R=basis.base_ring())

def gs_knh_fast_point(basis, degs, point, aux):
    """The function responsible for fixing a single point"""
    basis = gs_knh_fast_reduce_basis(basis, [point], aux)
    s, l = aux['s'], basis.dimensions()[0] - 1
    derivatives = [ (i,j) for i in range(0,s) for j in range(0,s-i) ] # order is important
    hasse = [ gs_Q_hasse_derivatives(basis[i], point, s, reduce=False) for i in range(l+1) ]
    PF = basis.base_ring(); F = PF.base_ring()
    T = identity_matrix(PF, l+1)
    for deriv in derivatives:
        evals = [ hasse[i][deriv] for i in range(l+1) ]
        # find minimal degree non-zero eval
        t = -1
        for i in range(l+1):
            if not evals[i].is_zero() and (t==-1 or degs[i] < degs[t]):
                t = i
        if t > -1:
            # fix T and hasse for other indices than hasse[t].
            coeffs = [ -evals[i]/evals[t] if  i != t and not evals[i].is_zero() else F.zero() for i in range(l+1) ]
            for i in range(l+1):
                if not coeffs[i].is_zero():
                    hasse[i] += hasse[t]*coeffs[i]
                    T.add_multiple_of_row(i, t, coeffs[i])
            T.set_row_to_multiple_of_row(t,t, x-point[0])
            # fix hasse[t] and degs[t]
            hasse[t] = zero_matrix(F, 1, s).stack(hasse[t].delete_rows([s-1]))
            degs[t] += 1
    return (T, degs)

def gs_knh_fast_tree(basis, degs, points, aux):
    """The D&C tree structure function of the fast KNH algorithm"""
    basis = gs_knh_fast_reduce_basis(basis, points, aux)
    if len(points)==1:
        return gs_knh_fast_point(basis, degs, points[0], aux)
    else:
        half = len(points)//2
        first  = points[:half]
        second = points[half:]
        (T, degs) = gs_knh_fast_tree(basis, degs, first, aux)
        basis = T * basis
        (T2, degs) = gs_knh_fast_tree(basis, degs, second, aux)
        return (T2 * T, degs)

def gs_construct_Q_knh_fast(points, tau, (s,l), wy):
    """Construct interpolation polynomial Q by the fast D&C variant of
    the Koetter--Nielsen--Hoeholdt algorithm"""
    F = points[0][0].parent()
    PF.<x> = F[]
    basis = identity_matrix(PF, l+1)
    degs = [ i*wy for i in range(l+1) ]
    aux = { 's': s, 'wy': wy, 'mods': gs_knh_fast_mod_tree(points, s, x) }
    (T, degs) = gs_knh_fast_tree(basis, degs, points, aux)
    # Resulting polys are basis*T = T
    PFy.<y> = PF[]
    # Construct Q as the element of the row with the lowest weighted degree
    codinglib.module.module_apply_weights(T, degs)
    best = codinglib.module.module_minimal_row(T)
    codinglib.module.module_remove_weights(T, degs)
    Qlist = T.row(best)
    return sum( y^i*PFy(Qlist[i]) for i in range(0,l+1) )







    
### Main Q construction if construction method does not matter
#Thresholds were determined by a set of tests for the paper
#  Fast Koetter-Nielsen-Hoeholdt interpolation in the G-S Algorithm
#     2014, Johan S. R. Nielsen, ACCT
#These tests were carried out on an Thinkpad T430s. I don't know how
#they vary over computers
#TODO: Instead of OR'ing on n and s, we could use a weighted average
#roughly corresponding to asymptotic complexity, say n^2*l^3*s
def gs_construct_Q(points, tau, (s,l), wy):
    if len(points) > 600 or s > 6:
        if len(points) > 2000 or s > 9:
            return gs_construct_Q_knh_fast(points, tau, (s,l), wy)
        else:
            return gs_construct_Q_knh(points, tau, (s,l), wy)
    else:
        return gs_construct_Q_lee_osullivan(points, tau, (s,l), wy)




#####################################################
### INSPECTION OF Q POLYNOMIALS
#####################################################

def gs_Q_hasse_derivatives(Q, point, s, reduce=True):
    """Compute all the Hasse derivatives of Q at point up to multiplicity s.
    If reduce==True then first compute Q mod (x-p[0])^s, and use this
    for Hasse derivative computation.
    """
    (a,b) = point
    F = a.parent()
    if hasattr(Q, "variables") and len(Q.variables()) == 2:
        # Q is probably in F[x,y]
        PF.<x> = F[]
        PFy.<y> = PF[]
        Q = PFy(Q)
    elif not hasattr(Q, "gen"):
        #assume Q is a list or vector of F[x]
        PF = Q[0].parent()
        x = PF.gen()
        PFy.<y> = PF[]
        Q = PFy(list(Q))
    else:
        PFy = Q.parent();
        PF = PFy.base_ring();
    F = PF.base_ring()
    x,y = PF.gen(), PFy.gen()
    if reduce:
        # Reduce Q to speed up computation
        G = (x-a)^s
        Qred = PFy([ Q[i] % G for i in range(Q.degree()+1) ]) 
        Qhatl = PFy(Qred(y=y+point[1], x=x+point[0])).list()
    else:
        Qhatl = PFy(Q(y=y+point[1], x=x+point[0])).list()
    def coeff(i,j):
        if i+j >= s or j>=len(Qhatl) or i > Qhatl[j].degree():
            return F.zero()
        else:
            return Qhatl[j].list()[i]
    return matrix(F, s, s, coeff)

def gs_Q_is_interpolation(Q, points, s):
    """Returns whether Q interpolates all the points with multiplicity s"""
    if Q.is_zero():
        return False
    for point in points:
        H = gs_Q_hasse_derivatives(Q, point, s)
        for dx in range(s):
            for dy in range(s-dx):
                if not H[dx,dy].is_zero():
                    return False
    return True

def gs_Q_weighted_degree(Q, wy):
    """Compute the weighted degree of a F[x][y] polynomial, where x weighs 1 and y weighs wy."""
    Ql = Q.list()
    return max( wy*i + (Ql[i].degree() if not Ql[i].is_zero() else -infinity) for i in range(len(Ql)))
    
#####################################################
### COMPLETE DECODER
#####################################################

class GRSDecoderGuruswamiSudan(DecoderBlockCodeAbstract):
    r"""
    Implementation of the Guruswami-Sudan decoding algorithm for GRS codes:
       Guruswami, V., and Sudan, M "Improved Decoding of Reed--Solomon Codes and
       Algebraic-Geometric Codes." IEEE Transactions on Information Theory 45,
       no. 6 (1999): 1757-67.

    INPUT:

    - ``C`` - The code to decode
    - ``tau=None`` - Optional, the sought decoding radius. Specify this or `params`
    - ``params=None`` - Optional, multiplicity `s` and list size `\ell` to use
      for decoding. Specify this or `tau`.
    - ``Qfinder=None`` - Optional, a function to use to find the Q polynomial.
      Must adhere to the calling convention of ``gs_construct_Q``, which is also
      used if None is given.

    FIELDS:
    
    - ``s`` - The multiplicity of the decoder
    - ``ell`` - The list size of the decoder
    - ``Qfinder`` - Function for finding Q(x,y). Defaults to `gs_construct_Q`
                     and must have same signature.
    - ``root_finder` - Function for finding y-roots of Q(x,y). Defaults to
                       `rootfinding.rootfind_bivariate` and must have same signature.
    """

    def __init__(self, C, tau=None, params=None, Qfinder=None, root_finder=None):
        self.C = C
        if tau:
            self.tau = tau
            self.s, self.ell = gs_minimal_list(C.n, C.k, tau)
        elif params:
            self.s = params[0]
            self.ell = params[1]
            (self.tau,_) = gs_decoding_radius(C.n, C.k, s=self.s, l=self.ell)
        else:
            raise Exception("Specify either tau or params")
        #TODO: Precompute stuff for the various Qfinders
        self.Qfinder = Qfinder if Qfinder else gs_construct_Q
        self.root_finder = root_finder if root_finder else codinglib.rootfinding.rootfind_bivariate

    def decoding_radius(self):
        return self.tau

    def is_list_decoder(self):
        return True

    def __str__(self):
        return "Guruswami-Sudan list-decoder up to %s for %s, (s,ell)=(%s,%s)" % (self.tau, self.C, self.s, self.ell)

    def decode_to_information(self, r, debug=0):
        if debug > 1:
            codinglib.codeTesting.time_ping()
        n,k,d,alphas,colmults = self.C.n, self.C.k, self.C.d, self.C.alphas, self.C.colmults
        ## SETUP INTERPOLATION PROBLEM
        wy = k-1
        points = [ (alphas[i], r[i]/colmults[i]) for i in range(0,len(alphas)) ]
        if debug > 1:
            codinglib.codeTesting.time_ping(mes="Setup")
        ## SOLVE INTERPOLATION
        Q = self.Qfinder(points, self.tau, (self.s,self.ell), wy)
        if debug > 1:
            codinglib.codeTesting.time_ping(mes="Solving interpolation")
        ## EXAMINE THE FACTORS AND CONVERT TO CODEWORDS
        factors = self.root_finder(Q, wy)
        if not factors:
            if debug > 0:
                print "No valid factors of interpolation polynomial"
            return None
        return [ poly2list(f, self.C.k) for f in factors ]
        #Note: there is no check that these are close to the received word.
