#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""2D key equations and solving them.

The name 2D Key Equations was introduced in [jsrn]. They are in Computer Algebra
sometimes known as Simultaneous Hermitian Pade approximations.

The notation here follows mainly that of [jsrn], but uses an improved handling
of weights described in [Dec].

- [jsrn]: Johan S. R. Nielsen, List decoding of algebraic codes, 2013 PhD thesis
          at Technical University of Denmark.
- [Dec]:  Johan S. R. Nielsen, Peter Beelen, Sub-quadratic Decoding of One-Point
          Hermitian Codes, in review for IEEE Transactions of Information Theory
"""
from codinglib.util import *
import codinglib.module

class KeyEquation2D():
    r"""Class for representing 2D key equations
    TODO: Currently only of Type 2
    TODO: Currently etas and ws are implemented in naive way. By making a
          shifted-polynomial class one can elegantly fix this

    INPUT:
    
    - ``sigma``  - Number of equations
    - ``rho``  - Number of $\Lambda$ variables
    - ``Gs``   - The moduli, a list of length $\sigma$ of $\mathbb F[x]$
    - ``Ss`` - The coefficients, a $\rho x \sigma$ matrix over $\mathbb F[x]$, or
      a list of $rho$ lists, each containing $\sigma$ polynomials.
    - ``nu``   - The x weight, must be positive integer.
    - ``etas`` - The $\Lambda$ weights, $\rho$ nonnegative integers
    - ``ws``   - The $\Omega$ weights, $\sigma$ nonnegative integers

    FIELDS:
    All of the above as well as

    - ``F``    - The base field
    - ``PF``   - The polynomial field
    - ``var``  - The variable of the polynomials
    """
    # Local fields
    # - ``_minimal_wmodule'': cached version of a minimised basis in weighted form

    def __init__(self, sigma, rho, Gs, Ss, etas, ws, nu=1, N=None, type=2):
        self.sigma = sigma
        self.rho = rho
        self.Gs = Gs
        if is_matrix(Ss):
            self.Ss = [ Ss.row(i) for i in range(0,Ss.nrows()) ]
        else:
            self.Ss = Ss
        self.nu = nu
        self.etas = etas
        self.ws = ws
        self.F = Gs[0].base_ring()
        self.PF = Gs[0].parent()
        self.var = self.PF.gen()
        self.type = type

        assert len(Gs) == sigma , "Specify exactly sigma Gs"
        assert len(self.Ss) == rho , "Specify exactly rho rows of Ss"
        assert all( len(Ssr)==sigma for Ssr in self.Ss ) , "Specify exactly sigma elements in each row of Ss"
        assert len(etas) == rho , "Specify exactly rho etas"
        assert len(ws) == sigma , "Specify exactly sigma ws"
        #TODO: Type parameter checks

    def get_module(self, weighted=False):
        """Calculate the $\mathbb F[x]$ matrix whose row space contains all the solutions"""
        def getMentry(r, c):
            if c < self.rho:
                return 1 if r==c else 0
            else:
                if r < self.rho:
                    return self.Ss[r][c-self.rho]
                else:
                    return self.Gs[r-self.rho] if r==c else 0
        M = Matrix(self.PF, self.rho + self.sigma, self.rho + self.sigma, getMentry)    
        if weighted:
            self.module_apply_weights(M)
            return M
        return M

    def module_weight_vector(self):
        """Return the weight vector for the module representation of the 2d key equation."""
        return list([ w/self.nu for w in self.etas + self.ws ])
        
    def module_apply_weights(self, M):
        """Assuming `M` is a basis of the module defined by this 2D key
        equation, apply inplace the column weights.

        The weighted basis is by embedding the weights and permuting the colums.
        The permutation for the columns is returned.
        See e.g. [Dec]"""
        return codinglib.module.module_apply_weights(M, self.module_weight_vector())

    def module_remove_weights(self, Mw):
        """Assuming `M_w` is a weighted basis of the module defined by this 2D
        key equation, return the unweighted basis"""
        codinglib.module.module_remove_weights(Mw, self.module_weight_vector())

    def module_minimised(self, weighted=False, algorithm='mulders'):
        """Return a minimal basis of the module defined by this 2d key
        equation. The result is cached and computed only once.

        Algorithm can be one of:
          - 'demand'  - Use the Demand--Driven algorithm
          - 'mulders' - Use the Mulders--Storjohann algorithm
        Note that even if algorithm is changed, the minimal basis will not be
        recomputed if it has already been computed.
        """
        if hasattr(self,"_minimal_wmodule"):
            if weighted:
                return copy(self._minimal_wmodule)
            else:
                M = copy(self._minimal_wmodule)
                self.module_remove_weights(M)
                return M
        if algorithm == 'demand':
            Lambdas = key2d_demand_driven(self)
            m = self.rho+self.sigma
            def get_row(i):
                # Complete the recorded Lambda solutions if non-zero. If zero, output the corresponding initial row.
                # In DD, row indexes equal leading position whenever i >= rho, which are the only possibility for
                # having zeroes on all Lambdas.
                try:
                    return self.solution_complete(Lambdas.row(i))
                except ValueError:
                    if i < self.rho:
                        raise Exception("Rows 0--(rho-1) should never have zeroes on all lambdas! Happened with row %s" % i)
                    return [ self.PF.zero() ] * i + [ self.Gs[i-self.rho] ] + [ self.PF.zero() ] * (self.rho+self.sigma-i-1)
            M = Matrix(self.PF, m, m, [ get_row(i) for i in range(self.rho+self.sigma) ])
            self.module_apply_weights(M)
            self._minimal_wmodule = M
        else: 
            #Use Mulders--Storjohann
            Mw = self.get_module(weighted=True)
            codinglib.module.module_weak_popov(Mw)
            self._minimal_wmodule = Mw
        return self.module_minimised(weighted=weighted)

    def solution_complete(self, lambdas):
        """Given the lambdas, compute the omegas which complete the solution,
        i.e. satisfy the congruences."""
        if all(lam.is_zero() for lam in lambdas):
            raise ValueError("All given Lambdas are zero")
        return vector(self.PF, self.rho+self.sigma,
                      [ lambdas[j] if j < self.rho
                        else sum( lambdas[h]*self.Ss[h][j-self.rho] for h in range(self.rho) ) % self.Gs[j-self.rho]
                        for j in range(self.rho+self.sigma) ])
        
    def solution_minimal(self, algorithm='mulders'):
        """Calculate a minimal solution by minimising the module. algorithm is like for ``self.module_minimised''."""
        mindeg, row = infinity, -1
        if algorithm=='mulders':
            M = self.module_minimised(weighted=False, algorithm=algorithm)
            weights = self.module_weight_vector()
            return M.row(codinglib.module.module_row_with_LP(M, list(range(self.rho)), weights=weights))
        else:
            raise NotImplementedError("demand-driven algorithm is currently out of order")
            #TODO: If using Demand-Driven, then it is wasteful to compute the entire basis if you want only the solution
                
    def solution_bound(self, weighted=False):
        r"""Return the bound for the maximal size of a minimal solution in the
        usual case. Proposition 2.39 of jsrn phd. If weighted=True, return the
        bound as the degree of the weighted smallest solution. Otherwise, return
        a vector of length $\rho+\sigma$ containing the unweighted degree
        constraints for each element.
        """
        if self.type == 1:
            raise NotImplementedError("Solution size bound for Type 1")
        else:
            #TODO: Corollary 2.40
            def lhs(d):
                return sum( pos(self.Gs[j].degree() - ceil((d-self.ws[j])/self.nu) ) for j in range(self.sigma) )
            def rhs(d):
                return sum(floor((d-eta)/self.nu) for eta in self.etas)  +self.rho - 1
            d = codinglib.util.find_minimal_satisfiable(lambda d: lhs(d) <= rhs(d), startn=1, contiguous=True)
            if weighted:
                return d
            else:
                return vector([ floor((d-self.etas[j])/self.nu) for j in range(self.rho) ]
                             +[  gilt((d-self.ws[j])/self.nu)   for j in range(self.sigma) ])

    def is_solution(self, v):
        """Test whether the vector $v$ is a solution or not. $v$ can be given as
        a complete solution, or only the lambdas."""
        if len(v) != self.rho and len(v) != self.rho+self.sigma:
            raise Exception("Invalid length for the purported solution")
        vt = self.solution_complete(v[:self.rho], weighted=False)
        if len(v) == self.rho + self.sigma:
            # Check the congruence equations are satisfied
            for j in range(0, self.sigma):
                if vt[j + self.rho] != v[j + self.rho]:
                    print "Congruence equations not satisfied: sigma_%s" % j
                    return False
        # Check degree constraints
        ldeg = max( self.nu * vt[i].degree() + self.etas[i] for i in range(0, self.rho) )
        odeg = max( self.nu * vt[j+self.rho].degree() + self.ws[j] for j in range(0, self.sigma) )
        if self.type == 1:
            return ldeg < self.N and odeg < self.N
        else:
            return ldeg > odeg




### Operating functions        
def key2d_previous(ke, theta, h):
    r"""For the set of weights specified by the 2D key equation, return the legal
    degree and leading position which attain the maximal possible value less
    than that of theta and h.
    Note that if $theta \notequiv w_h \mod h$, then the result will probably not
    be legal."""
    if ke.nu == 1:
        if h > 0:
            return (theta, h-1)
        else:
            return (theta-1, ke.rho+ke.sigma-1)
    else:
        if hasattr(ke, "_previous_set"):
            pset = ke._previous_set
        else:
            # Build previous_set
            ws = ke.etas + ke.ws
            m  = len(ws)
            # sort cols prim. by high weight mod nu, secondly by right-most first
            hs = list(range(m))
            hs.sort(key = lambda h: -( m*(ws[h] % ke.nu) + h ) )
            pset = dict()
            for j in range(0, m):
                next = j+1 if j < m-1 else 0
                dtheta = (ws[hs[j]] - ws[hs[next]]) % ke.nu
                pset[hs[j]] = (hs[next], dtheta)
            ke._previous_set = pset
                
        (h2, dtheta) = pset[h]
        return (theta-dtheta, h2)


def key2d_demand_driven(ke):
    """Solve the 2D key equation using the demand--driven algorithm (of jsrn
    phd). Returns the first $rho$ columns of a basis in weak Popov form."""
    #TODO: Fast coefficient getter when Gs[h] is sparse
    #TODO: nu is represented naively
    raise NotImplementedError("Demand-driven solving of 2d key equations is currently out of order")
    Mw = ke.get_module(weighted=True)
    x = ke.var
    Lambda = Mw[:, 0:ke.rho]
    thetas = []
    alphas = []
    def swap_rows(i, h):
        Lambda.swap_rows(i, h)
        thetas[i], thetas[h] = thetas[h], thetas[i]
        alphas[i], alphas[h] = alphas[h], alphas[i]
    # Create coefficient-getter functions, for preprocessing sparsity of the Gs etc.
    def _get_coefficient(h):
        def direct_coefficient(p, theta):
            if theta > p.degree():
                return 0
            else:
                return p[theta]
        if h < ke.rho:
            return lambda i, thetai: direct_coefficient(Lambda[i, h], thetai)
        else:
            Gw = Mw[h, h]
            return lambda i, thetai: \
                direct_coefficient( sum( Lambda[i, j].shift(-ke.etas[j])*Mw[j, h] for j in range(ke.rho) )
                                         % Gw, thetai)
    get_coefficient = [ _get_coefficient(h) for h in range(ke.rho+ke.sigma) ]

    # The actual algorithm, close to line-by-line equivalent to pseudo-code
    W = set()
    for i in range(ke.rho + ke.sigma):
        lp = LP(Mw.row(i))
        thetas.append(Mw[i, lp].degree())
        alphas.append(Mw[i, lp].leading_coefficient())
        if lp != i:
            W.add((i, lp))
    while W:
        (i, h) = W.pop()
        if thetas[i] < thetas[h]:
            swap_rows(i, h)
        Lambda.add_multiple_of_row(i, h, -alphas[i]/alphas[h]*x^(thetas[i]-thetas[h]))
        rep = True
        while rep:
            (thetas[i], h) = key2d_previous(ke, thetas[i], h)
            alphas[i] = get_coefficient[h](i, thetas[i])
            rep = alphas[i].is_zero()
        if i != h:
            j = None
            for pair in W:
                if pair[0] == h:
                    j = pair[1]
                    break
            if j:
                swap_rows(i, h)
                W.remove((h, j))
                W.add((i, j))
            else:
                W.add((i, h))
    codinglib.module.module_remove_weights(Lambda, ke.etas)
    #TODO: _remove_exponent has been removed
    return _remove_exponent(Lambda, ke.nu)
