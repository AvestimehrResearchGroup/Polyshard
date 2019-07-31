#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Class for representing GRS codes with corresponding utility functions"""
from codinglib.util import *
from codinglib.code import BlockCodeAbstract
import codinglib.codeTesting

class GRS(BlockCodeAbstract):
    r"""
    Class for evaluation-style GRS codes
    
    INPUT:

    - ``F`` - The field of the code
    - ``n`` - The length of the code
    - ``k`` - The dimension of the code
    - ``alphas`` - The evaluation points. Must be n different elements of F
    - ``colmults=None`` - Column multipliers. Must be n long. Assumed all 1.

    FIELDS:

    - ``F`` - Field
    - ``n`` - Code length
    - ``k`` - Code dimension
    - ``d`` - Minimum distance
    - ``alphas`` - Evaluation points
    - ``colmults`` - Column multipliers
    """
    def __init__(self, F, n, k, alphas, colmults=None):
        assert n==len(alphas)  , "The no. of evaluation points alphas must be n"
        assert len(alphas)==len(set(alphas)) , "the evaluation points should all be distinct"
        assert all( a in F for a in alphas ) , "the evaluation points should all be in F"
        assert k>0 and k<=n , "The dimension k should be between 1 and n"
        if not colmults is None:
            assert len(colmults) == n , "The no. of column multipliers must be n"
            assert all(c in F for c in colmults) , "The column multipliers should all be in F"

        self.F = F
        self.n = n
        self.k = k
        self.d = n-k+1
        self.alphas = alphas
        self.colmults = colmults if not colmults is None else [F.one()]*n

    def __str__(self):
        return "[%s,%s,%s] GRS code over %s" % (self.n,self.k,self.d,self.F)

    def __repr__(self):
        return self.__str__()

    def _latex_(self):
        return r"[%s,%s,%s]{\rm\ GRS\ code\ over\ } %s" % (self.n,self.k,self.d,latex(self.F))

    def __eq__(self,other):
        """Checks for total equivalence of two GRS codes, i.e. completely same parameters."""
        return isinstance(other, GRS) \
           and self.F == other.F \
           and self.alphas == other.alphas \
           and self.k == other.k \
           and self.colmults == other.colmults
        
    @cached_method
    def generator_matrix(self):
        """Return the generator matrix for the GRS evaluation code."""
        return matrix(self.F, self.k, self.n,
                lambda i,j: self.colmults[j]*self.alphas[j]^i \
                    if self.alphas[j] != 0 else (0 if i!=0 else 1))

    def syndrome(self, r, direct=True):
        """Calculate the syndrome polynomial from the received word (as a
        list). If ``direct`` is False, it is the syndrome fitting together with
        the indirect error locator, otherwise the direct."""
        ws = self.parity_colmults()
        if not direct:
            ws = self.parity_colmults()
            ls = []
            ap = [ self.F.one() ] * self.n
            for i in range(0, self.d-1):
                ls.append(sum([ r[j]*ws[j]*ap[j] for j in range(0, self.n) ]))
                if i < self.d-2:
                    for j in range(self.n):
                        ap[j] *= self.alphas[j] 
            return ls
        else:
            P = x.parent()
            return [ sum( r[j]*ws[j]*self.alphas[j]^(self.d-i-2)
                        for j in range(0,self.n) )
                        for i in range(0,self.d-1) ]

    def error_locator(self, errVec, x, direct=True):
        """Calculate the error locator in the variable x given the error
        positions. If direct, calculate the direct error locator."""
        errPos = codinglib.codeTesting.error_pos(errVec)
        if not direct:
            return prod([ 1-x*self.alphas[i] for i in errPos ])
        else:
            return prod([ x-self.alphas[i] for i in errPos ])
        
    def error_evaluator(self, errVec, x, direct=True):
        """Calculate the error evaluator in the variable x given the error
        positions and error values. If direct, calculate the direct error
        locator."""
        errPos = codinglib.codeTesting.error_pos(errVec)
        if not direct:
            return sum( [ errVec[i]*prod([ 1-x*self.alphas[j]
                            for j in errPos if j!=i ])
                            for i in errPos ] )
        else:
            ws = self.parity_colmults()
            return sum( (-e[i]*self.alphas[i]^(self.d-1)*ws[i])
                        *prod( x-self.alphas[j] for j in errPos if j!=i )
                    for i in errPos )

    def error_evaluator_factor(self, i, errPos, direct=True):
        """The error evaluator evaluated at ``alpha_i'' returns the
        error value multiplied by a certain factor. This function returns that
        factor"""
        if not direct:
            raise NotImplementedError
        else:
            ws = self.parity_colmults()
            a = self.alphas
            return -a[i]^(self.d-1)*ws[i] \
                        *prod( a[i]-a[h] for h in errPos if h != i )

    def reencode_received_word(self, r, positions):
        """Reencode the given received word at the given ``k'' positions: this
        means adding a codeword to ``r'' such that the result is 0 at those
        positions"""
        if positions != range(self.k):
            raise NotImplementedError("Currently, only reencoding at first k positions is supported")
        PF = PolynomialRing(self.F, 'x')
        f_r = PF.lagrange_polynomial(zip(self.alphas[:self.k], r[:self.k]))
        return vector([ self.F.zero() ]*self.k + [ r[i] - f_r(self.alphas[i]) for i in range(self.k,self.n) ])

    @cached_method
    def multipliers_product(self):
        r"""The product of the j'th column multiplier and the j'th parity
        check column multiplier. This is the inverse of the products of the
        $\alpha_j-\alpha_i$ for $i\neq j$."""
        a = self.alphas
        return [ 1/prod([ a[i] - a[h] for h in range(0, len(a)) if h != i ])
                    for i in range(0,len(a)) ]

    @cached_method
    def parity_colmults(self):
        """The column multipliers of the parity check matrix, that is, the
        column multipliers of the dual code. 
        After first calculation, they are cached, so just call this."""
        n = self.n
        cs = self.colmults
        etas = self.multipliers_product()
        return [ etas[i]/cs[i] for i in range(0,n) ]

    def low_weight_codeword(self, weight=None):
        """Return a random low-weight codeword from self. If weight is not
        given, return a minimum-weight codeword (i.e. weight $n-k+1$).
        Otherwise, return one of weight if possible.
        NOTE: The codeword returned is only truly random amongst legally
        returnable words if weight = n-k+1"""
        # Find the word by constructing an evaluation polynomial which has
        # exactly n-weight zeroes amongst the alphas
        if not weight:
            weight = self.n-self.k+1
        if weight < self.n-self.k+1 or weight > n:
            raise Exception("weight must be between n-k+1 and n, both inclusive")
        from random import shuffle
        zeroes = [ a for a in self.alphas ] # copy list alphas
        shuffle(zeroes)
        zeroes = zeroes[:self.n-weight]
        F = self.F
        P.<x> = F[]
        p = P.one() * prod( x-zero for zero in zeroes )
        a = F.zero()
        while not a:
            a = F.random_element()
        #TODO: Multiply p with a random poly with no zeroes in alphas\zeroes.
        info_poly = a*p
        info = list(info_poly) + [F.zero()]*(self.k - p.degree()-1)
        return self.encode(info)

    def true_minimum_distance(self):
        """The true minimum distance of a GRS code is always $n-k+1$"""
        return self.n-self.k+1

    def true_dimension(self):
        """The true dimension of a GRS code is always $n-k+1$"""
        return self.k


class RSFrequency(GRS):
    r"""
    Conveniency class for RS codes defined using the frequency domain language.
    These are GRS codes where the length divides q-1 and the evaluation points
    are all powers of an n'th root of unity in decreasing order.

    For such codes, the indirect syndrome calculation is exceedingly simple, as
    it is just the evaluation of the received word, seen as a polynomial, for
    all n-k higher powers of the n'th root of unity. For this reason, these
    codes default to using the indirect error locator.
    
    INPUT:

    - ``F`` - The field of the code
    - ``n`` - The length of the code
    - ``k`` - The dimension of the code
    - ``alphas`` - The evaluation points. Must be n different elements of F

    FIELDS: As for a GRS code.
    """

    def __init__(self, F, n, k):
        assert( (F.cardinality()-1) % n == 0 )
        alpha = F.multiplicative_generator()^((F.cardinality()-1)/n)
        alphas = [ alpha^(-i) for i in range(0, n) ]
        #Cdual = GRS(F, n, k, alphas)
        #GRS.__init__(self, F, n, k, alphas, colmults=Cdual.parity_colmults())
        colmults = [ F(1/n) for i in range(0, n) ]
        GRS.__init__(self, F, n, k, alphas, colmults=colmults)

    def syndrome(self, r, direct=False):
        return GRS.syndrome(self, r, direct)

    def error_locator(self, errVec, x, direct=False):
        """Calculate the error locator in the variable x given the error
        positions. If direct, calculate the direct error locator."""
        return GRS.error_locator(self, errVec, x, direct)
        
    def error_evaluator(self, errVec, x):
        """Calculate the error evaluator in the variable x given the error
        positions and error values. If direct, calculate the direct error
        locator."""
        return GRS.error_evaluator(self, errVec, x, direct=False)

    def error_evaluator_factor(self, i, errVec, direct=False):
        """The error evaluator evaluated at ``alpha_i'' returns the
        error value multiplied by a certain factor. This function returns that
        factor"""
        return GRS.error_evaluator_factor(self, i, errVec, direct)
