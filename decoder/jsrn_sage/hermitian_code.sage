#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Class for representing and decoding Hermitian codes.

The codes follow the definition of [jsrn], which is a slightly larger class of
codes than in [Dec].
The decoding algorithms are described in [Dec].

- [jsrn]: Johan S. R. Nielsen, List decoding of algebraic codes, 2013 PhD thesis
          at Technical University of Denmark.
- [Dec]:  Johan S. R. Nielsen, Peter Beelen, Sub-quadratic Decoding of One-Point
          Hermitian Codes, in review for IEEE Transactions of Information Theory
"""

#TODO: Method for point groups

from codinglib.util import *
from codinglib.code import *
from codinglib.key2d import KeyEquation2D
from codinglib.module import *
from codinglib.rootfinding import rootfind_modulo

from codinglib.codeTesting import clock_call

class HermitianCode(BlockCodeAbstract):
    r"""
    Class for evaluation-style one-point Hermitian codes, with evaluation points chosen in groups.

    The form used for the Hermitian curve is `y^q + y = x^{q+1}` over a field
    `F` of cardinality `q^2` for some `q`.
    
    Choose integers `n, m` and affine, rational places `P_1,\ldots,P_n` on the
    Hermitian curve, with the `P_i` chosen in groups, such that for any
    `F`-element, the `q` places with this same `x`-coordinate (written in
    coordinate style) are either all chosen or none of them are. Let `P_\infty`
    be the rational place at infinity. The code is then

    `C = { f(P_1),\ldots,f(P_n) \mid f \in L(mP_\infty) }`

    Currently supports only `m > 2g-2` where `g = 1/2 q (q-1)` is the genus of the curve.
    
    INPUT:

    - `q` - The square-root of the cardinality of the field of the code
    - `m` - The pole order of evaluation functions.
    - `alphas` - Optional: the evaluation point groups, members of `F_{q^2}`.
      Each alpha gives q evaluation points. If not specified, than all elements
      `F_{q^2}` are assumed.

    FIELDS:

    - `q` - The square-root of the cardinality of the field of the code
    - `F` - Field of the code; has cardinality `q^2`
    - `n` - Code length
    - `k` - Code dimension
    - `d` - Designed minimum distance = `n-m`
    - `alphas` - Evaluation point groups
    - `g` - the genus of the corresponding Hermitian curve
    - `Ps`- evaluation points, members of `F_{q^2}^2` and on the Hermitian
      curve.
    - `eval_powers` - pairs (i,j) such that `order(x^i*y^j) <= m`, i.e.
      monomials used for evaluation. They are used for evaluation in the order
      listed here
    - `Ya` - the ring of functions with only poles at P_infinity, i.e.
      polynomials in x and y. It's degree function is minus the pole order at
      infinity, but only valid for polynomials of y-degree at most q. This is
      represented simply as a bivariate polynomial ring in x and y with a
      corresponding monomial order, so it does *not* automatically consider the
      curve equation between x and y.
    """

    def __init__(self, q, m, alphas=None):
        self.F = GF(q^2,'a')
        if alphas==None:
            alphas = self.F.list()
        else:
            assert all(alpha in self.F for alpha in alphas), "All alphas must belong to F_{q^2}"

        Ps = []
        for alpha in alphas:
            aq = alpha^(q+1)
            count = 0
            for beta in self.F: 
                if beta^q + beta == aq:
                    Ps.append((alpha, beta))
                    count += 1
            assert ( count == q )
        self.Ps = Ps

        self.g = (q*(q-1))//2
        assert m > 2*self.g - 2, "Only codes with m > 2g - 2 = %s are supported." % (2*self.g-2)
        self.q = q
        self.m = m
        self.alphas = alphas
        self.k = m - self.g + 1
        self.n = len(Ps)
        assert self.n == q * len(alphas), "Something is fishy with the evaluation groups."
        self.d = self.n-m

        self.Ya = PolynomialRing(self.F,'x,y',order=TermOrder('wdeglex',(q,q+1)))

        monoms = []
        for i in range(0, m//q+1):
            xd = i * q
            for j in range(0, min(q-1, (m-xd)//(q+1)) +1 ):
                monoms.append((i,j))
        self.monoms = monoms
        assert(len(monoms) == self.k)

    def __str__(self):
        return "[%s,%s,>=%s] Hermitian code over %s" \
                    % (self.n,self.k,self.d,self.F)

    def __repr__(self):
        return self.__str__()

    def _latex_(self):
        return r"[%s,%s,%s]{\rm\ Hermitian\ code\ over\ } %s" % (self.n,self.k,self.d,latex(self.F))
        
    @cached_method
    def generator_matrix(self):
        """Return the generator matrix for the Hermitian code, evaluation style."""
        def cell(row, col):
            P = self.Ps[col]
            (i,j) = self.monoms[row]
            return power_discr(P[0], i) * power_discr(P[1], j)
        return matrix(self.F, self.k, self.n, cell)

    def information_to_function(self, mes):
        """Return the function in F[x,y] which corresponds to the information
        vector ``mes``, as returned by ``self.unencode(c)`` for a codeword
        ``c``"""
        assert len(mes) == self.k , "The message vector does not have the length of the dimension of the code"
        return self.Ya({ self.monoms[i]: mes[i] for i in range(len(self.monoms)) })

    def function_to_information(self, f):
        """Return the information vector in F which corresponds to the function f in F[x,y]"""
        assert f.degree() <= self.m , "f has too high pole order at infinity"
        x,y = self.Ya.gens()
        return vector( f.monomial_coefficient(x^i*y^j) for (i,j) in self.monoms )

    def true_dimension(self):
        """The true dimension of a Hermitian code is always $m - g + 1$"""
        return self.k

    def is_error_locator(self, Lambda, err_vec):
        """Returns whether the given function ``Lambda`` in ``Ya`` locates the errors of ``err_vec``"""
        if Lambda.is_zero():
            return all( ei.is_zero() for ei in err_vec)
        for i in range(self.n):
            if not(err_vec[i].is_zero()) and not(Lambda(x=self.Ps[i][0], y=self.Ps[i][1]).is_zero()):
                return False
        return True




class HermitianDecoderAbstract(DecoderBlockCodeAbstract):
    r"""Abstract base class for decoding algorithms for Hermitian codes,
    advertising some specialised helper functions.

    Several helper functions are for conversion between the three
    representations of elemens of `Ya`: Standard basis (F[x,y] element with
    y-degree at most q-1); Shat basis (F[x,y] element with x-degree at most q);
    and power series in x.

    To achieve the complexities that are discussed in [Dec] for the conversion
    functions, we rely on the Sage in-built support for sparse polynomials (uni-
    and multi-variate). Outside the conversion functions themselves, the
    polynomials are always represented densely. Inside the conversion functions,
    and in particular in some of the cached/precomputed values, the sparse
    polynomials are saved as dictionaries.
    """

    def __init__(self, C):
        self.C = C
        #Dense polynomial ring in x
        Px = PolynomialRing(C.F, 'x')
        self.Px = Px
        #Sparse polynomial ring in x
        sPx = PolynomialRing(C.F, 'x', sparse=True)
        self.sPx = sPx
        #Sparse polynomial ring in x,y
        sPxy = PolynomialRing(C.F, 'x,y', sparse=True)
        self.sPxy = sPxy

    #@clock_call
    def precompute(self):
        self._lagrange_function_precompute()

    @cached_method
    def _monomial_to_sparse_power_series_d(self, i, j, N):
        r"""Like _monomial_to_sparse_power_series but returns an element the
        sparse polynomial ring self.sPx"""
        x = self.sPx.gen()
        if i>0:
            zero_x_poly = self._monomial_to_sparse_power_series_d(0,j,N)
            return x^i*zero_x_poly

        # i == 0 now.
        if j==0:
            return self.sPx.one()
        if j==1:
            # y = sum_{b=0}^\infty (-1)^b x^{(q+1)q^b}
            q = self.C.q
            one = self.sPx.one()
            maxb = floor(log(N/(q+1), q))
            return self.sPx({ (q+1)*q^b: (-one)^b for b in range(maxb+1) })
        else:
            # Do things recursively to avoid ever getting too large powers
            half = j//2
            t1 = self._monomial_to_sparse_power_series_d(0,half,N)
            t2 = self._monomial_to_sparse_power_series_d(0,j-half,N)
            return (t1*t2).truncate(N)

    @cached_method
    def _monomial_to_sparse_power_series(self, i, j, N):
        r"""Compute the power series in x of x^i*y^j up to precision N, i.e. x^N
        is not included (sparse).

        This is returned as a dictionary (as given by self.sPx.dict())"""
        return self._monomial_to_sparse_power_series_d(i,j,N).dict()

    #@clock_call
    def _std_basis_to_power_series(self, p, N):
        r"""Convert `p \in \Ya` represented as element in F[x,y] to power series of precision up to N (dense)"""
        q = self.C.q
        pl = [self.C.F.zero()]*N
        for (coeff, mon) in p:
            i,j = mon.degrees()
            gij = self._monomial_to_sparse_power_series(i, j, N)
            for t in gij:
                if t < N:
                    pl[t] += coeff * gij[t]
        return self.Px(pl)

    @cached_method
    def _y_power_to_sparse_std_d(self, j):
        r"""Like _y_power_to_sparse_std but result represented in the sparse `F[x,y]` ring."""
        q = self.C.q
        b,a = j.quo_rem(q) 
        (sx,sy) = self.sPxy.gens() 
        if a+b <= q-1:
            if a>0:
                return sy^a * self._y_power_to_sparse_std_d(j-a)
            else:
                return (sx^(q+1) - sy)^b
        else:
            binom = self._y_power_to_sparse_std_d(b*q)
            p1 = binom.truncate(sy, q-a)
            p2 = binom//sy^(q-a)
            return sy^a*p1 + (sx^(q+1) - sy)*p2

    @cached_method
    def _y_power_to_sparse_std(self, j):
        r"""Convert `y^j` to an element in `F[x,y]` with `y`-degree less than `q`.
        This is represented as a dictionary (as given by self.sPxy.dict()).
        Note that it is assumed that `j < q^2`."""
        return self._y_power_to_sparse_std_d(j).dict()

    #@clock_call
    def _shat_basis_to_std(self, p):
        r"""Convert `p \in \Ya` represented in `F[x,y]` using the Shat basis (low
        x-degree) into one using the S basis (low y-degree) (both dense).
        Note that it is assumed that y-degree of `p` is less than q^2.
        """
        q = self.C.q
        (x,y) = p.parent().gens() 
        m = p.degree()
        assert p.degrees()[1] < q^2 , "y-degree of p is too large"
        sp = self.sPxy(p.truncate(y, q)).dict() #Take everything in S cap Shat, represent as dict
        pred = p//y^q
        for (coeff, mon) in pred: #Go through the remaining monomials
            (i, jred) = mon.degrees()
            # sp += coeff*gj    as dicts
            gj = self._y_power_to_sparse_std(jred + q) 
            for gmon in gj:
                mon = gmon.eadd_p(i,0)
                if mon in sp:
                    sp[mon] += coeff*gj[gmon]
                else:
                    sp[mon] = coeff*gj[gmon]
        #assert sp.degrees()[1] < q , "y-degree was not reduced to less than q!"
        return self.C.Ya(sp) # Convert back into dense

    #@clock_call
    def _power_series_to_shat_basis(self, p, N):
        r"""Convert `p \in \Ya` represented as a power series in x into an
        element of F[x,y] in the Shat basis (dense)"""
        # Go through x^i*y^j in order of increasing (0,0) order of vanishing
        q = self.C.q
        mons = dict()
        pl = p.list()
        pl = pl + [p.parent().zero()]*(N-len(pl))
        for h in range(N):
            if not pl[h].is_zero():
                j,i = Integer(h).quo_rem(q+1)
                gij = self._monomial_to_sparse_power_series(i,j,N)
                coeff = pl[h] / gij[h]
                for t in gij:
                    if t<N:
                        pl[t] -= coeff * gij[t]
                mons[(i,j)] = coeff
        return self.C.Ya(mons)

    def _power_series_to_std_basis(self, p, N):
        r"""Convert `p \in \Ya` represented as a power series in x into an
        element of F[x,y] in the standard basis (dense)"""
        return self._shat_basis_to_std(self._power_series_to_shat_basis(p, N))

    def _xlist_to_Ya(self, Pls):
        r"""Take a a list of `y`-coefficients, i.e. polys in `x`, and return a polynomial in `Ya`"""
        y = self.C.Ya.gens()[1]
        return sum(y^i*self.C.Ya(Pls[i]) for i in range(len(Pls)))

    def _Ya_to_xlist(self, p):
        r"""Take a polynomial in `F[x,y]` and return a list of its `q` `y`-coefficients, i.e. polys in `x`"""
        y = self.C.Ya.gens()[1]
        if p.degree(y) >= self.C.q:
            raise Exception("Ya-element not properly reduced")
        return vector( self.Px( p.coefficient({y: j})) for j in range(self.C.q))

    def _xlist_multiply_by_y(self, a):
        r"""Multiply `a \in \F[x,y] represented as an xlist by `y`"""
        x = self.Px.gen()
        return [a[-1]*x^(self.C.q+1) , a[0] - a[-1]] + [ a[i-1] for i in range(2, self.C.q) ]

    def _multiply_and_reduce(self, a, b):
        r"""For `a, b \in F[x,y]` both with `y`-degree less than `q`, compute
        `ab` with representation having `y`-degree less than `q`.

        This implementation is supposed to be slightly faster than using the
        _shat_to_std_basis, which would actually work."""
        q = self.C.q
        x,y = self.C.Ya.gens()
        ab = a*b
        if ab.degrees()[1] < q:
            return ab
        ab_tr = ab.truncate(y, q) #low degrees
        ab_top = ab//y^q           #high degrees
        return ab_tr + (x^(q+1) - y)*ab_top

    @cached_method
    def _lagrange_function_alpha_unit(self, alphai):
        (x,y) = self.C.Ya.gens()
        return prod( (x-alpha)/(alphai-alpha) for alpha in self.C.alphas if not alphai==alpha)

    @cached_method
    def _lagrange_function_beta_unit(self, (alphai, betai)):
        (x,y) = self.C.Ya.gens()
        B = []
        for (alpha,beta) in self.C.Ps:
            if alpha==alphai:
                B.append(beta)
        return prod( (y-beta)/(betai-beta) for beta in B if not betai==beta)

    def _lagrange_function_slow(self, r):
        r"""Given a vector `r`, compute an `R in \Ya` such that `R(P_i) = r_i` for `i=1,\ldots,n`, iterative strategy"""
        C = self.C
        (x,y) = self.C.Ya.gens()
        P = self.C.Ya.zero()
        for i in range(self.C.n):
            P += r[i] * self._lagrange_function_alpha_unit(C.Ps[i][0]) * self._lagrange_function_beta_unit(C.Ps[i])
        return P

    def _lagrange_function_precompute(self):
        """Precompute various stuff for computing Ya Lagrange interpolation polys."""
        #Construct the point map: a map `d: Fqq -> (Fqq -> int)` such that
        #d[alpha][beta] = i whenever (alpha,beta) is the i'th point of the code
        C = self.C
        pm = dict()
        def add(alpha, beta, i):
            if not alpha in pm:
                pm[alpha] = dict()
            Bs = pm[alpha]
            Bs[beta] = i
        for i in range(C.n):
            P = C.Ps[i]
            add(P[0], P[1], i)
        self._lagrange_function_point_map = pm
        A = list(pm.keys())
        self._lagrange_function_A = A

        # compute multiplier on all r[i] for which C.Ps[i]= (alpha, *)
        self._lagrange_function_mults = { alpha: 1/prod( alpha - alphap for alphap in A if alpha != alphap ) for alpha in A } 

        # make sure all xprods have been cached
        def x_rec(Abegin,Aend):
            if Abegin < Aend:
                mid = (Aend-Abegin)//2 + Abegin
                x_rec(Abegin, mid) ; x_rec(mid+1, Aend)
                self._lagrange_function_xprod(Abegin,mid) ; self._lagrange_function_xprod(mid+1,Aend)

    @cached_method
    def _lagrange_function_xprod(self, Abegin, Aend):
        A = self._lagrange_function_A
        (x,_) = self.C.Ya.gens()
        return prod( x-A[i] for i in range(Abegin, Aend+1) )

    def _lagrange_function_fast(self, r):
        r"""Given a vector `r`, compute an `R in \Ya` such that `R(P_i) = r_i` for `i=1,\ldots,n`, D&C strategy"""
        C = self.C
        F = C.F
        if not hasattr(self,'_lagrange_function_A'):
            self._lagrange_function_precompute()
        pm = self._lagrange_function_point_map
        A = self._lagrange_function_A
        mults = self._lagrange_function_mults
        Ry = PolynomialRing(F, 'y')
        def lagr_rec(Abegin,Aend):
            """Interpolate through all points with alpha-coordinate from Abegin,..,Aend"""
            if Abegin < Aend:
                mid = (Aend-Abegin)//2 + Abegin
                p1 = lagr_rec(Abegin, mid)
                p2 = lagr_rec(mid+1, Aend)
                return p1 * self._lagrange_function_xprod(mid+1,Aend) + p2 * self._lagrange_function_xprod(Abegin,mid)
            else:
                # Abegin == Aend
                alpha = A[Abegin]
                B = pm[alpha]
                g = mults[alpha] * Ry.lagrange_polynomial(list( (beta, r[B[beta]]) for beta in B ))
                return C.Ya(g)
        return lagr_rec(0, len(A)-1)

    def _lagrange_function(self, r):
        r"""Given a vector `r`, compute an `R in \Ya` such that `R(P_i) = r_i` for `i=1,\ldots,n`."""
        P = self._lagrange_function_fast(r)
        assert P.degrees()[1] < self.C.q , "y-degree of Lagrange should be at most q-1"
        assert P.degree() < self.C.n + 2*self.C.g , "Order of Lagrange is %s but should be less than n + 2g = %s" % (P.degree(), self.C.n + 2*self.C.g )
        return P



class HermitianDecoderPower(HermitianDecoderAbstract):

    def __init__(self, C, ell):
        HermitianDecoderAbstract.__init__(self, C)
        self.ell = ell
        assert ell * C.m < C.n , "We require ell * m < n"
        Px = PolynomialRing(C.F, 'x')

    def precompute(self):
        HermitianDecoderAbstract.precompute(self)
        x = self.Px.gen()
        self.G = self.Px( prod( x-alpha for alpha in C.alphas ) )
        
    def is_list_decoder(self):
        return False

    def decoding_radius(self):
        """The probabilistic bound on the decoding radius, i.e. the number of errors usually correctable"""
        #TODO: Is this what we want?
        ell = self.ell
        return min( floor( ell/(ell+1) * self.C.n - ell*self.C.m/2 - ell/(ell+1) ) , self.C.n - self.C.m - self.C.g )

    def decoding_radius_assured(self):
        """The number of errors which is guaranteed to be correctable"""
        return (self.C.d-1-self.C.g)//2

    def __str__(self):
        return "Power decoder with ell=%s for %s" % (self.ell, self.C)

    def decode_to_information(self, r):
        if not hasattr(self, 'G'):
            self.precompute()
        Lambda, Lambda_f = self._find_solution_pair(r)
        f = self._attempt_division(Lambda, Lambda_f)
        if f:
            return self.C.function_to_information(f)
        else:
            raise DecodingFailedError("Division of Lamdba*f and Lambda did not work")

    def _find_solution_pair(self, r):
        C = self.C
        Rs = self._power_lagranges(r)
        KE = self._build_key_equation(Rs)
        s = self._lattice_reduction(KE)
        return self._xlist_to_Ya(list(s[:C.q])) , self._xlist_to_Ya(list(s[C.q:2*C.q]))

    #@clock_call
    def _lattice_reduction(self, KE):
        # Mw =KE.get_module(weighted=True)
        # module_weak_popov(Mw)
        # KE._minimal_wmodule = Mw
        return KE.solution_minimal(algorithm='mulders')

    #@clock_call
    def _power_lagranges(self, r):
        Rs = list( self._lagrange_function(vector( power_discr(ri, t) for ri in r )) for t in range(1,self.ell+1) )
        assert(all( Rs[t+1](x=self.C.Ps[i][0], y=self.C.Ps[i][1]) == power_discr(r[i], t) for i in range(self.C.n))
               for t in range(self.ell)  )
        return Rs

    #@clock_call
    def _build_key_equation(self, Rs):
        C = self.C
        x = self.Px.gen()
        rho, sigma = C.q, C.q*self.ell
        Rls = [ self._Ya_to_xlist(Ri) for Ri in Rs ]
        def cell(ii,jj):
            (i, j) = ii+1, jj+1
            t = floor((j-1)/C.q) + 1
            h = (j-1) % C.q + 1
            if h>1:
                if i < h:
                    return Rls[t-1][h-i]
                elif i == h:
                    return Rls[t-1][0] - Rls[t-1][C.q-1]
                else:
                    return ( -Rls[t-1][C.q-1-(i-h)] + x^(C.q+1)*Rls[t-1][C.q-(i-h)] ) % self.G
            else:
                if i==1:
                    return Rls[t-1][0]
                else:
                    return x^(C.q+1)*Rls[t-1][C.q - i+1] % self.G
        Ss = Matrix(self.Px, rho, sigma, cell)

        etas = [ (i-1)*(C.q+1) + self.ell*C.m + 1  for i in range(1, C.q+1) ]
        ws = []
        for j in range(1, sigma+1):
            t = floor((j-1)/C.q) + 1
            h = (j-1) % C.q + 1
            ws.append( (h-1)*(C.q+1) + (self.ell-t)*C.m )
        return KeyEquation2D(sigma, rho, [self.G]*sigma, Ss, nu=C.q, etas=etas, ws=ws)

    #@clock_call
    def _attempt_division(self, Lambda, Lambda_f):
        N = self.C.q^3 
        Ndelta = 2*N #TODO: Optimise
        pLambda = self._std_basis_to_power_series(Lambda, Ndelta)
        pLambda_f = self._std_basis_to_power_series(Lambda_f, Ndelta)
        x = pLambda.parent().gen()
        if pLambda.is_zero():
            raise Exception("Lambda was zero")
        delta = 0
        while pLambda[delta].is_zero():
            delta += 1
        pLambda = (pLambda//x^delta).truncate(N)
        pLambda_f = (pLambda_f//x^delta).truncate(N)

        (g,pLambda_inv,b) = pLambda.xgcd(x^N)
        assert g == pLambda.parent().one() , "GCD or delta not correctly computed\n: 1 != %s" % g

        pf = (pLambda_f * pLambda_inv).truncate(N)
        f = self._power_series_to_std_basis(pf, N)
        if f.degree() <= self.C.m:
            return f
        else:
            return None
        # assert f.degree() <= self.C.m , "Degree of the found function after division is too high: %s" % f.degree()



class HermitianDecoderGuruswamiSudan(HermitianDecoderAbstract):
    __MES_IMPOSSIBLE_PARAMS = "Impossible parameters for the Guruswami-Sudan algorithm"

    def __init__(self, C, tau=None, params=None):
        HermitianDecoderAbstract.__init__(self, C)
        if tau:
            self.tau = tau
            self.s, self.ell = HermitianDecoderGuruswamiSudan.get_parameters_satisfiable(C, tau)
        else:
            self.s, self.ell = params
            self.tau = HermitianDecoderGuruswamiSudan.decoding_radius_from_parameters(C, params)
        Px = PolynomialRing(C.F, 'x')
        self.Px = Px
        self.Z = PolynomialRing(C.Ya, 'z')
        self.PxZ = PolynomialRing(Px, 'z')

    def precompute(self):
        HermitianDecoderAbstract.precompute(self)
        x = self.Px.gen()
        self.G = self.Px( prod( x-alpha for alpha in self.C.alphas ) ) # = x^(q^2) - x when all points were chosen
        self.Gs = [ 1, self.G ] + [self.Px.zero()]*(self.ell-1)
        for t in range(2, self.ell+1):
            self.Gs[t] = self.G * self.Gs[t-1]

    ### Static class methods for parameter calculation

    @staticmethod
    def _coefficient_count_Ya(C, T):
        """Return the no. of available coefficients in Ya-elements of degH < T."""
        if not hasattr(C, "_gap_counts"):
            gaps = set()
            gap_counts = [0 for i in range(2*C.g+1)]
            for i in range(1,2*C.g+1):
                if ((i - C.q) < 0 or (i-C.q) in gaps) and ((i-C.q-1) < 0 or (i-C.q-1) in gaps):
                    gaps.add(i)
                    gap_counts[i] = gap_counts[i-1]+1
                else:
                    gap_counts[i] = gap_counts[i-1]
            C._gap_counts = gap_counts
        if T > 2*C.g:
            return T - C.g
        elif T < 0:
            return 0
        else:
            return T - C._gap_counts[T]

    @staticmethod
    def are_parameters_satisfiable(C, tau, (s, ell)):
        """Return whether the given parameters are satisfiable for Guruswami-Sudan decoding the given Hermitian code C."""
        if s <= 0 or ell <= 0:
            return False
        T = s*(C.n - tau)
        coeffs = sum( HermitianDecoderGuruswamiSudan._coefficient_count_Ya(C, T - j*C.m) for j in range(0, ell+1))
        restrictions = C.n*s*(s+1)//2
        return coeffs > restrictions

    @staticmethod
    def get_parameters_satisfiable(C, tau):
        """Compute parameters (s,ell) for Guruswami-Sudan decoding the Hermitian
        code C, given the sought decoding radius.

        This implementation is similar to gs_minimal_list in gs.sage, in that it
        finds the minimal s parameters using binary search (as opposed to a closed
        formula in e.g. gs_params).
        There is a small complication since the number of coefficients available
        for a given s,ell is not a nice, smooth expression, so after an initial
        rough, binary-search estimate, we proceed by linear attempts
        """ 
        def try_s(s):
            # Safeguard against infinite loop
            if s > C.n^2:
                raise ValueError("The requested decoding radius seems to be impossible")
            (min_ell,max_ell) = \
                    solve2deg_int(C.m , C.m + 2*C.g - 2*s*(C.n-tau) , C.n*s*(s+1) + 2*C.g - 2*s*(C.n-tau) )
            if max_ell >= 0 and max_ell >= min_ell:
                return max(1, min_ell) 
            else:
                return None
        s = find_minimal_satisfiable(try_s, 1, contiguous=True)
        ell = try_s(s)
        # This pair s, ell works, but see if we can reduce s a bit
        while True:
            if   HermitianDecoderGuruswamiSudan.are_parameters_satisfiable(C, tau, (s-1,ell-1)):
                s,ell = s-1,ell-1
            elif HermitianDecoderGuruswamiSudan.are_parameters_satisfiable(C, tau, (s-1,ell)):
                s = s-1
            elif HermitianDecoderGuruswamiSudan.are_parameters_satisfiable(C, tau, (s,ell-1)):
                ell = ell-1
            else:
                break

        assert HermitianDecoderGuruswamiSudan.are_parameters_satisfiable(C, tau, (s,ell)) ,\
            "Logical error in parameter calculation: this should never happen."
        return (s,ell)

    @staticmethod
    def decoding_radius_from_parameters(C, (s, ell)):
        """Compute the decoding radius of Guruswami-Sudan decoding the Hermitian code C, given the parameters s and ell."""
        #TODO: This does not take into account the improved decoding radius analysis of Peter Beelen.
        # In particular, I believe the -C.g/s term should be -C.g/s/2.
        if s <= 0 or ell <= 0:
            return 0
        #The following is a pretty good estimate, but is not always correct
        tau = gilt( C.n*(1 - (s+1)/2/(ell+1)) - ell*C.m/2/s - C.g/s )
        # Attempt to increase by exact calculation
        while HermitianDecoderGuruswamiSudan.are_parameters_satisfiable(C, tau+1, (s,ell)):
            tau += 1
        return tau


    

    ### Normal object methods

    def is_list_decoder(self):
        return True

    def decoding_radius(self):
        """The assured bound on the decoding radius, i.e. the number of errors we can always correct."""
        return HermitianDecoderGuruswamiSudan.decoding_radius_from_parameters(self.C, (self.s, self.ell))

    #TODO: Add the "probably decode up to" radius?

    def __str__(self):
        return "Guruswami-Sudan decoder with (s,ell)=(%s,%s) for %s" % (self.s, self.ell, self.C)

    def decode_to_information(self, r):
        if not hasattr(self, 'G'):
            self.precompute()
        C = self.C
        M = self._build_matrix(r)
        Q = self._construct_Q_from_matrix(M)
        assert self._Z_degree(Q) < self.s*(C.n - self.tau) \
            , "Something's wrong: the found Q has too large degree: %s" % (self._Z_degree(Q))
        ls = self._root_find(Q)
        infos = []
        for f in ls:
            if f.degree() <= C.m:
                #TODO: Possibly check for distance of encoded word
                infos.append(C.function_to_information(f))
        if infos:
            return infos
        else:
            return None

    ### HELPER FUNCTIONS

    def _Z_multiply_by_Ya(self, H, p):
        """Compute `pH` where `H \in Ya[z]` and `p \in Ya`, all `Ya`-elements
        represented in the standard basis."""
        return self.Z([self._multiply_and_reduce(Hi, p) for Hi in H])

    def _Z_to_Ya_list(self, H):
        r"""Convert `H \in \Ya[z] to a list of \Ya elements of length `\ell+1`"""
        return H.list() + [self.C.Ya.zero()]*(self.ell+1 - (H.degree()+1))

    def _xlist_to_Ya_list(self, xlist):
        r"""For a list of `F[x]` elements of length `q(\ell+1)`, representing an
        element in `\Ya[z]`, bundle its elements into a list of `\ell+1` `Ya` elements."""
        C, q = self.C, self.C.q
        ya_list = [ C.Ya.zero() ]*(self.ell+1)
        for t in range(self.ell+1):
            ya_list[t] = self._xlist_to_Ya(xlist[t*q : (t+1)*q])
        return ya_list

    def _Ya_list_to_Z(self, ya_list):
        r"""The inverse of 'self._Z_to_Ya_list'."""
        return self.Z(ya_list)

    def _Z_degree(self, Q):
        r"""Return `\deg_{H,m}(Q)`, using the notation of [Dec]."""
        Qlist = Q.list()
        return max( Qlist[t].degree() + t*self.C.m for t in range(len(Qlist)) )

    def _Z_to_PxZ(self, Q, k):
        r"""Convert `Q \in \Ya[z]` into `\F[[x]][z]` of power series precision `k`."""
        return self.PxZ(list(self._std_basis_to_power_series(Qi, k) for Qi in Q.list()))



    ### MAIN STEPS OF THE ALGORITHM

    #@clock_call
    def _build_matrix(self, r):
        C = self.C
        # Calculate the Ht, first as Ya[z]-elems, then as lists of xlists (i.e. ell+1 lists of q F[x]-elems each)
        z = self.Z.gen()
        R = self._lagrange_function(r)
        zRs = [ self.Z.one(), z-R ]
        for t in range(2,self.s+1):
            zRs.append( z*zRs[t-1] - self._Z_multiply_by_Ya(zRs[t-1], R) )
        Hs = [ self._Z_multiply_by_Ya(zRs[t], self.Gs[self.s-t]) for t in range(self.s+1) ] \
             + [ z^(t-self.s)*zRs[self.s] for t in range(self.s+1, self.ell+1) ]
        Hs = [ self._Z_to_Ya_list(Ht) for Ht in Hs ]
        Hs = [ [ self._Ya_to_xlist(Ht_ya) for Ht_ya in Ht] for Ht in Hs ]
        # For each Ht, add the q y-shifts of it as a row to the matrix
        (x,y) = C.Ya.gens()
        rows = []
        for Hi in Hs:
            h = Hi # this is ell+1 lists of q F[x] elems each
            for t in range(0, C.q):
                rows.append(list(flatten_once(h)))
                if t < C.q-1:
                    for i in range(self.ell+1):
                        h[i] = self._xlist_multiply_by_y(h[i])
        return matrix(self.Px, C.q*(self.ell+1), C.q*(self.ell+1), rows)

    #@clock_call
    def _construct_Q_from_matrix(self, M):
        C, q = self.C, self.C.q
        weights = list(flatten_once([ (t*C.m + i*(q+1))/q for i in range(q) ] for t in range(self.ell+1)))
        module_mulders_storjohann(M, weights=weights)
        rowi = module_minimal_row(M, weights)
        return self._Ya_list_to_Z(self._xlist_to_Ya_list(M.row(rowi)))

    #@clock_call
    def _root_find(self, Q):
        C = self.C
        k = self._Z_degree(Q) + 1
        Qpow = self._Z_to_PxZ(Q, k)
        # Get roots as a list of truncated power series
        pow_roots = rootfind_modulo(Qpow, precision=k)
        assert all( d_h > C.m for (_,d_h) in pow_roots) , "NotImplemented: Handling roots that do not have sufficient precision"
        return [ self._power_series_to_std_basis(root, C.m+1) for (root, d_h) in pow_roots ]
