#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Class for representing Goppa codes with corresponding utility
functions"""

from codinglib.util import *
from codinglib.code import *

class GoppaBinary(BlockCodeAbstract):
    r"""
    Class for evaluation-style binary Goppa codes.

    INPUT:

    -  ``bigF`` - The extension field of the code

    -  ``n`` - The length of the code

    -  ``g`` - The generator polynomial in ``bigF[x]``

    -  ``alphas`` - The evaluation points. A list of different elements from
                    ``bigF``.

    FIELDS:

    -  ``F`` - F2
    - ``bigF`` - extension field
    - ``g`` - Goppa polynomial
    - ``m`` - size of the extension ``[F:F2]`` (``log2`` of the cardinality)
    - ``n`` - code length
    - ``k`` - lower bound on code dimension. At least 0.
    - ``d`` - lower bound on minimum distance
    - ``alphas`` - evaluation points
    """

    def __init__(self, bigF, n, g, alphas):
        if not alphas:
            alphas = F.list()[1:] #Take all non-zero elements
        self.F = GF(2)
        self.bigF = bigF
        self.m = log(bigF.cardinality(), 2)
        self.n = n
        self.g = g
        self.k = max(n - self.m*g.degree(), 0)
        self.d = 2*g.degree()+1
        self.alphas = alphas
        #TODO: Sanity check the parameters

    def __str__(self):
        return "[%s,>=%s,>=%s] binary Goppa code over %s" \
                    % (self.n,self.k,self.d,self.F)

    def __repr__(self):
        return self.__str__()

    def _latex_(self):
        return r"$[%s,>=%s,>=%s]${\rm\ binary\ Goppa\ code\ over\ } %s" \
                    % (self.n,self.k,self.d,latex(self.bigF))

    def __eq__(self,other):
        """Checks for total equivalence of two Goppa codes, i.e. completely same parameters."""
        return isinstance(other, GoppaBinary) \
           and self.bigF == other.bigF \
           and self.alphas == other.alphas \
           and self.g == other.g
        

    @cached_method
    def parity_check_matrix_big_field(self):
        """Return a parity check matrix for this code. Uses the formula (12) of
        MacWilliams & Sloane, p. 340."""
        r = self.g.degree()
        X = matrix(self.bigF, r, self.n,
                lambda i,j: self.alphas[j]^i if i != 0 else 1)
        Y = matrix(self.bigF, self.n, self.n,
                lambda i,j: 1/self.g(self.alphas[i]) if i==j else 0)
        H = X*Y
        return H

    @cached_method
    def parity_check_matrix(self):
        """Return a parity check matrix for this code. Uses the formula (12) of
        MacWilliams & Sloane, p. 340."""
        H = self.parity_check_matrix_big_field()
        return matrix_subfield(H)

    def syndrome(self, r):
        """Calculate the Goppa syndrome polynomial from word r (as a list)."""
        (x,) = self.g.variables()
        FQuot = self.g.parent().quotient_ring(self.g)
        Squot = sum( FQuot(r[i])/FQuot(x-self.alphas[i]) for i in range(0,self.n) )
        return Squot.lift()

    def error_locator(self, errPos, x):
        """Calculate the error locator in the variable x given the error
        positions. Always the direct error locator."""
        return prod([ x-self.alphas[i] for i in errPos ])
        
    def error_evaluator(self, errPos, x):
        """Calculate the error evaluator in the variable x given the error
        positions. Always the direct error evaluator."""
        return sum( prod(x-self.alphas[j] for j in errPos if j!=i)
                        for i in errPos )
