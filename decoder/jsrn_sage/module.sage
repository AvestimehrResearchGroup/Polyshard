#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

r"""Functionality for free F[x]-modules, in particular for calculating
bases in weak Popov form. The weak Popov form is slightly stronger
than row reduced form, and it is also a type of reduced Groebner bases
over some simple module monomial orderings.

- [Dec]:  Johan S. R. Nielsen, Peter Beelen, Sub-quadratic Decoding of One-Point
          Hermitian Codes, in review for IEEE Transactions of Information Theory
"""

from codinglib.util import *

def normalise_weights(weights):
    """Return the normalisation of the weight vector, i.e. add/subtract t from
    all weights, such that all weights are non-negative and at least one is 0"""
    m = min(weights)
    return [ w-m for w in weights ]

################################################################################
# Module minimisation / Lattice basis reduction algorithms
################################################################################

def module_row_reduction(M, i, j, pos):
    """Perform a row reduction with row j on row i, reducing position
    pos. If M[i,pos] has lower degree than M[j,pos], nothing is changed.
    Returns the multiple of row j used."""
    pow = M[i,pos].degree() - M[j,pos].degree()
    if pow < 0:
        return None
    coeff = -M[i, pos].leading_coefficient() / M[j, pos].leading_coefficient()
    x = M.base_ring().gen()
    multiple = x^pow*coeff
    M.add_multiple_of_row(i, j, multiple)
    return x^pow*coeff

def module_perform_row_reduction(M, last=None):
    """Perform some row transformation to M in-place if possible: find some conflict and
    reduce it. A conflict is two rows having the same leading position. If last
    is given, it should be a possible conflict which will be used if it is an
    actual conflict."""
    conflict=None
    if last:
        if LP(M.row(last[0])) == LP(M.row(last[1])):
            conflict = last
    if not conflict:
        def findConflict(M, weights=None):
            """If M is a matrix over F[x], find two rows which have the same
            leading position"""
            lpos = [ LP(M.row(i), weights) for i in range(0, M.nrows()) ]
            posmap = dict()
            for i in range(0,M.nrows()):
                lposi = lpos[i]
                if lposi>-1:
                    if lposi in posmap:
                        return (i, posmap[lposi])
                    else:
                        posmap[lposi] = i
            return None
        conflict = findConflict(M)
    if conflict:
        (i,j) = conflict
        lpos = LP(M.row(i))
        # Make sure i is the greater column
        if M[i,lpos].degree() < M[j,lpos].degree():
            (i,j) = (j,i)
        return (i,j, module_row_reduction(M, i, j, lpos))
    return None

def module_mulders_storjohann(M, weights=None, debug=0):
    """Reduce $M$ to weak Popov form using the Mulders--Storjohann
    algorithm (Mulders, T., and A. Storjohann. "On Lattice Reduction
    for Polynomial Matrices." Journal of Symbolic Computation 35, no.
    4 (2003): 377-401.).

    Handles column weights with only a slight possible penalty. The weights can
    be fractions, but cannot be floating points.
    
    If debug is True, then print some info."""
    if weights and len(weights) != M.ncols():
        raise ValueError("The number of weights must equal the number of columns")
    # initialise conflicts list and LP map
    LP_to_row = dict( (i,[]) for i in range(M.ncols()))
    conflicts = []
    for i in range(M.nrows()):
        lp = LP(M.row(i), weights=weights)
        ls = LP_to_row[lp]
        ls.append(i)
        if len(ls) > 1:
            conflicts.append(lp)
    iters = 0
    # while there is a conflict, do a row reduction
    while conflicts:
        lp = conflicts.pop()
        ls = LP_to_row[lp]
        i, j = ls.pop(), ls.pop()
        if M[i,lp].degree() < M[j, lp].degree():
            j,i = i,j

        module_row_reduction(M, i, j, lp)
        ls.append(j)
        lp_new = LP(M.row(i), weights=weights)
        if lp_new > -1:
            ls_new = LP_to_row[lp_new]
            ls_new.append(i)
            if len(ls_new) > 1:
                conflicts.append(lp_new)
        iters += 1
    return iters



################################################################################
# Higher-level reduction algorithms
################################################################################

def module_weak_popov(M, weights=None, debug=0):
    """Compute a (possibly column-weighted) weak Popov form of $M$.
    The weights can be any non-negative fractions."""
    return module_mulders_storjohann(M, weights=weights, debug=debug)

def module_popov(M, weights=None):
    """Reduce M to the unique (row-wise) Popov form of M and return the
    permutation matrix needed for achieving Popov form from weak Popov form.

    A matrix is in (row-wise) Popov form if it is in weak Popov form, and,
    furthermore, that
    1) every leading position is monic.
    2) The leading positions are in ascending order, counted by rows in ascending order
    3) All entries ``M[h, i]`` have degree less than ``M[i,LP(i)]``.
    """
    module_weak_popov(M, weights=weights)

    # print module_all_degrees(M, weights)
    # print 

    nrows = M.nrows()
    #Make monic
    lpos = [ LP(M.row(i), weights) for i in range(0, nrows) ]
    for i in range(0,nrows):
        M.set_row_to_multiple_of_row(i,i,1/M[i,lpos[i]].leading_coefficient())
    #Sort the rows by leading position
    #Note: presumably to satisfy Permutation, matrix.permute_rows is 1-indexed instead of 0-indexed!
    if weights is None:
        weights = [0]*M.ncols()
    ipos = [ (i+1, lpos[i]) for i in range(M.nrows()) ]
    ipos.sort(key = lambda (_,lp): lp)
    pi = Permutation([ i for i,_ in ipos ])
    M.permute_rows(pi)

    # print module_all_degrees(M, weights)
    
    #Perform reductions with earlier rows on later ones, in the right order
    # Starting with the highest-degree rows
    lpos = [ LP(M.row(i), weights) for i in range(0, nrows) ]
    idegs = [ (i, M.ncols()*(M[i, lpos[i]].degree() + weights[lpos[i]]) + lpos[i]) for i in range(M.nrows()) ]
    idegs.sort(key = lambda (_,d): -d)
    degree_order = [ i for i,_ in idegs ]
    changed = True
    while changed: # fingers crossed. I didn't prove this doesn't give infinite loop
        changed = False
        for i in degree_order:
            for j in range(0, nrows):
                if j != i:
                    ideg = M[i, lpos[i]].degree() + weights[lpos[i]]
                    # Cancel all other positions in that column having high degree
                    while M[j,lpos[i]].degree() + weights[lpos[i]] >= ideg:
                        module_row_reduction(M,j,i,lpos[i])
                        changed = True
            #             print "\t", i, j, lpos[i]
            # print i
            # print module_all_degrees(M, weights)
    # Make sure everything is all right
    assert module_is_popov(M, weights) , "Internal error: Matrix is still not Popov"
    return pi


def module_order_basis(M, d, weights=None):
    """Compute the left order basis of `M` up to degree `d`.
    
    This is a matrix `F` that is a row reduced basis for all vectors `v` such
    that `vM` is zero modulo `x^d`.
    """
    Px = M.base_ring()
    x = Px.gen()
    xd = x^d
    n, m = M.dimensions()
    if weights:
        if len(weights) != n:
            raise ValueError("Weights must be a list/vector of length M.nrows()")
        weights = normalise_weights(weights)
        max_weight = max(weights)
        weightsB = weights + [d+1+max_weight]*m
    else:
        weightsB = [0]*n + [d+1]*m
    def cellB(i,j):
        if j<n:
            return Px.one() if i==j else Px.zero()
        if i<n:
            return M[i, j-n]
        return xd if i==j else Px.zero()
    B = matrix(Px, n+m, n+m, cellB)
    module_weak_popov(B, weights=weightsB)
    null_rows = []
    for i in range(B.nrows()):
        row = B.row(i)
        if row[n:].is_zero():
            null_rows.append(row[:n])
    if len(null_rows) != n:
        raise Exception("Error: The order basis computation did not produce the expected number of null rows")
    return matrix(Px, n,n, null_rows)
     

################################################################################
# Query algorithms 
################################################################################

def module_is_row_reduced(M, weights=None):
    """Return whether the matrix is row reduced.

    Row reduced means that the leading position matrix of M has the same rank as M.

    Note: This function is quite slow since it uses Sage functions to compute
    the rank of `M`. In many cases, it would probably be faster to reduce `M` to
    weak Popov form and assert that
    """
    if module_is_weak_popov(M, weights):
        return True
    return module_LP_matrix(M, weights).rank() == module_rank(M, weights)
 
def module_is_weak_popov(M, weights=None):
    """Return whether the matrix is in weak Popov form"""
    seen = set()
    for i in range(M.nrows()):
        lp = LP(M.row(i), weights=weights)
        if lp in seen:
            return False
        seen.add(lp)
    return True

def module_is_popov(M, weights=None):
    """Returns true if M is in Popov form, otherwise throws an Assertion
    error.

    A matrix is in (row-wise) Popov form if it is in weak Popov form, and,
    furthermore, that
    1) every leading position is monic.
    2) The leading positions are in ascending order, counted by rows in ascending order
    3) All entries ``M[h, i]`` have degree less than ``M[i,LP(i)]``.
    """
    nrows, ncols = M.nrows(), M.ncols()
    lpos = [ LP(M.row(i), weights) for i in range(0, nrows) ]
    for i in range(0, nrows):
        # leading position should be monic
        if not M[i,lpos[i]].leading_coefficient().is_one() :
            # "The %s'th row does not have monic leading position" % i
            return False
        # leading position should be increasing
        if i > 0:
            if lpos[i] < lpos[i-1]:
                return False
                # "The rows are not properly sorted (at row %s)" % i
        # the leading position should have greatest degree of that column
        for j in range(0, nrows):
            if j!=i:
                if M[i, lpos[i]].degree() < M[j, lpos[i]].degree():
                    # "The %s'th row's leading position (%s) is not dominant in that column" % (i, lpos[i])
                    return False
    return True

def module_LP_matrix(M, weights=None):
    """Return the leading position matrix of M.

    If M is over F[x], then the leading position matrix L is over F and has the
    same dimensions as M. L[i,j] is then the coefficient to `x^{d_i}` in M[i,j],
    where `d_i` is the (weighted) degree of the `i`'th row of M.
    """
    P = M.base_ring()
    F = P.base_ring()
    ldegs = [ module_degree_of_row(M, i, weights) for i in range(M.nrows()) ]
    def cell(i,j):
        deg = ldegs[i]
        if weights:
            deg -= weights[j]
        return 0 if deg > M[i,j].degree() else M[i,j][deg]
    L = matrix(F, M.nrows(), M.ncols(), cell)
    return L

def module_rank(M, weights=None):
    """Compute the rank of a module M.

    We do this by reducing it to weak Popov form and reading off the number of
    non-zero rows.
    """
    if module_is_weak_popov(M, weights):
        Mc = M
    else:
        Mc = copy(M)
        module_weak_popov(Mc)
    return len( [ i for i in range(Mc.nrows()) if not Mc.row(i).is_zero() ])

def module_degree(M, weights=None):
    """Return the maximal degree occurring in M."""
    return max( module_degree_of_row(M, i, weights) for i in range(M.nrows()) for j in range(M.ncols()))

def module_all_degrees(M, weights=None):
    """Return a matrix containing the degrees of all the entries of M"""
    if not weights:
        return matrix(SR, M.nrows(), M.ncols(), lambda i,j: M[i,j].degree() if not M[i,j].is_zero() else -infinity)
    else:
        return matrix(SR, M.nrows(), M.ncols(), lambda i,j: M[i,j].degree() + weights[j] if not M[i,j].is_zero() else -infinity)

def module_degree_of_row(M, i, weights=None):
    """Return the degree of the i'th row of M, i.e. the greatest degree among
    its elements."""
    lp = LP(M.row(i), weights)
    if lp < 0:
        return None
    lt = M.row(i)[lp]
    return None if lt.is_zero() else (lt.degree() + weights[lp] if weights else lt.degree())

def module_row_degrees(M, weights=None):
    """Return the row-degrees of $M$ as a vector over ZZ."""
    return vector(ZZ, [module_degree_of_row(M, i, weights=weights) for i in range(0, M.nrows())])

def module_orthogonality_defect(M):
    """Return the orthogonality defect of the square matrix $M$ over a
    polynomial ring. Requires the calculation of the determinant of $M$."""
    return module_row_degree(M) - M.determinant().degree()

def module_construct_vector(Mp, v):
    """If the polynomial vector v lies in Mp, where Mp is a module by
    rows, and is in weak Popov form, then return the F[x]-linear combination of
    the rows of Mp which create v. This is essentially the multivariate division
    algorithm."""
    rowMap = dict([ (LP(Mp.row(i)) , i) for i in range(0,Mp.nrows()) ])
    combs = [Mp.base_ring().zero() for t in range(0, Mp.nrows())]
    x = Mp.base_ring().gen()
    while not v.is_zero():
        pos = LP(v)
        redPoly = Mp[rowMap[pos], pos]
        pow = v[pos].degree() - redPoly.degree()
        if pow < 0:
            return None
        coeff = v[pos].leading_coefficient() / redPoly.leading_coefficient()
        multiple = coeff * x^pow
        v = v - multiple*Mp.row(rowMap[pos])
        combs[rowMap[pos]] += multiple
    return combs

def module_contains(Mp, v):
    """Test whether the polynomial vector v lies in Mp, where Mp is a module by
    rows, and is in weak Popov form."""
    return module_construct_vector(Mp, v) != None

def module_subset_module(M1, M2):
    """Test whether M2 is a subset of M1 as modules by rows, by transforming a
    copy of M1 to weak popov form and asserting that the bases of M2 are in this."""
    M1c = copy(M1)
    module_weak_popov(M1c)
    for i in range(0, M2.nrows()):
        if not module_contains(M1c, M2.row(i)):
            print "Row %s failed" % i
            return False
    return True

def module_fractional_weight_permutation(weights):
    """Return the permutation which can be used for embedding the semantics of
    fractional weights into the module minimisation.
    A permutation is returned of the integers from 1 to ``len(numerators)``.
    See e.g. [Dec].
    """
    n = len(weights)
    denominator = lcm(list(f.denominator() for f in weights))
    numerators = [ f.numerator() * denominator/f.denominator() for f in weights ]
    residues = [ num % denominator for num in numerators ]
    res_map = dict()
    for i in range(n):
        if residues[i] in res_map:
            res_map[residues[i]].append(i+1)
        else:
            res_map[residues[i]] = [i+1]
    res_uniq = sorted(res_map.keys())
    return Permutation(list(flatten_once([ res_map[res] for res in res_uniq ])))
        
def module_apply_weights(M, weights):
    """Apply column weights inplace to the matrix `M`. If weights are all integers, then
    `M` is multiplied on the $n$th column with `x^{weights[n]}`.
    If weights are fractions, then `M` is appropriately column permuted and
    multiplied with the `x^t` where `t = int(weights[n])`. Afterwards, the
    permutation is returned if needed for reference. See e.g. [Dec]
    """
    if all( w.is_integer() for w in weights):
        for j in range(M.ncols()):
            M.set_col_to_multiple_of_col(j,j, x^weights[j])
    else:
        perm = module_fractional_weight_permutation(weights)
        for j in range(M.ncols()):
            M.set_col_to_multiple_of_col(j,j, x^floor(weights[j]))
        M.permute_columns(perm)
        return perm

def module_remove_weights(M, weights):
    """Remove the weights as introduced by ``module_apply_weights``."""
    #M = M.change_ring(QQ)
    if all( w.is_integer() for w in weights):
        for i in range(M.nrows()):
            for j in range(M.ncols()):
                M[i,j] = M[i,j].shift(-weights[j])
    else:
        perm = module_fractional_weight_permutation(weights)
        pinv = perm.inverse()
        M.permute_columns(pinv)
        module_remove_weights(M, [ floor(wj) for wj in weights ])

def module_row_with_LP(M, lp, weights=None):
    """If LP is a column index, return any row with this column as leading
    position. If lp is a list of column indices, return the row with minimal
    degree (possibly weighted) and has leading position in this set.
    """
    lp = set(lp if hasattr(lp, '__iter__') else [lp])
    (best,besti) = infinity, -1
    for i in range(M.nrows()):
        r = M.row(i)
        lpi = LP(r, weights=weights)
        if lpi in lp:
            d = module_degree_of_row(M, i, weights=weights) #-1 if r[lpi].is_zero() else (r[lpi].degree() + weights[lpi] else  r[lpi].degree())
            if (not d is None) and (d < best or best is None):
                (best,besti) = d, i
    return besti

def module_minimal_row(M, weights=None):
    """Return the (least) index of the (non-zero) row with minimal degree"""
    (best,besti) = module_degree_of_row(M, 0, weights), 0
    for i in range(1, M.nrows()):
        d = module_degree_of_row(M, i, weights)
        if (not d is None) and (d < best or best is None):
            (best,besti) = d, i
    return besti
            
def module_maximal_row(M, weights=None):
    """Return the (least) index of the (non-zero) row with maximal degree"""
    (best,besti) = module_degree_of_row(M, 0, weights), 0
    for i in range(1, M.nrows()):
        d = module_degree_of_row(M, i, weights)
        if (not d is None) and (d < best or best is None):
            (best,besti) = d, i
    return besti
