#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"Cyclic linear codes."

from codinglib.util import *
from codinglib.code import BlockCodeAbstract

class CyclicCode(BlockCodeAbstract):
    r"""
    Class for polynomial style cyclic codes.

    INPUT:

    - ``F`` - Field
    - ``n`` - code length
    - ``g=None`` - generator polynomial
    - ``D=None`` - Defining set of the code (might be partially defined)

    Either g or D must be supplied at construction. If D is supplied, the
    actual defining set used is the smallest set possible such that the
    resulting generator matrix fits in F. This is done by taking the cyclotomic
    coset closure of D.

    FIELDS:

    - ``F``      - field
    - ``n``      - code length
    - ``k``      - code dimension
    - ``g``      - generator polynomial
    - ``Fsplit`` - the splitting field of g
    - ``s``      - The log_q of the cardinality of Fsplit
    - ``alpha``  - The nth root of unity in Fsplit used for the defining set 
    - ``D``      - the complete defining set of the code, i.e. all powers of alpha which are zeroes of g.
    - ``cosets`` - powers i of F.gen such that D = union(M_i for i in cosets)
      where M_i is the cyclotomic coset corresponding to the i'th power of
      F.gen.
    
    METHODS:

    In excess of the usual Code methods, CyclicCode supports the following
    methods.

    - ``coset_poly(coset)`` -  The coset polynomial corresponding to the cycl.
      coset. Is in F[x]. Note that g = prod(coset_poly(d) for d in D)
    - ``ell()`` - The smallest integer such that n | (F.cardinality()^ell - 1)
    - ``minimum_distance_BCH()`` - The BCH bound of this code
    - ``minimum_distance_HT()`` - The Hartmann-Tzeng bound of this code
    """
    
    
    def __init__(self, F, n, g=None, D=None):
        self.F = F
        q = F.cardinality()
        self.n = n
        if not g and not D:
            raise Exception("Either g or D must be specified")

        s = find_modulo_power(n, -1, q)
        if not s:
            raise Exception("No such s such that n | q^s - 1")
        self.s = s
        try:
            Fsplit.<beta> = GF(q^s, 'beta', repr='log')
        except TypeError: #If repr is unsupported
            Fsplit.<beta> = GF(q^s, 'beta')
        self.Fsplit = Fsplit
        alpha = beta^((q^s-1)/n) #primitive nth root of unity in Fsplit
        self.alpha = alpha
        if g:
            raise NotImplementedError("Defining cyclic code using gen poly")
            self.g = g
            (x,) = g.variables()
            assert ((x^n-1) % g).is_zero() , "Generator does not divide x^n-1"


            roots = Fsplit(g).roots()
            assert all(f.degree()==1 for (f,m) in roots)

            #TODO
            roots = set(self.D) # copy
            cosets = []
            while roots:
                r = roots.pop()
                cosets.append(r)
                roots = roots.difference(cyclotomic_coset(r))
        else: # D was defined
            cosets = []
            pows = set()
            for r in D:
                if not r in pows:
                    cosets.append(r)
                    pows = pows.union(cyclotomic_coset(n, r, q))
            self.cosets = cosets
            self.D = pows
            PFsplit.<xx> = Fsplit['x']
            PF = F['x']
            self.g =  PF( prod( xx-alpha^p for p in pows ) )
            self.k = n-self.g.degree()
        
    @cached_method
    def generator_matrix(self):
        """Return the generator matrix for this cyclic code. This is n-r-1 rows
        with the coefficients of g(x) shifted one place for each row."""
        r = self.g.degree()
        coeffs = self.g.list()
        G = matrix(self.F, self.k, self.n,
                lambda i,j: coeffs[j-i] if i<=j and j-i<=r else 0)
        return G

    #TODO: Should input/output list and not polynomial to conform with standard
    def encode(self, info):
        """Calculate the codeword corresponding to the info information word
        polynomial or list"""
        if isinstance(info, list):
            (x,) = self.g.variables()
            infop = list2poly(info, x)
        else:
            infop = info
        return (info*self.g) % (x^self.n-1)

    def __str__(self):
        #TODO: improve
        return "[%s,%s,?] cyclic code over %s" % (self.n, self.k, self.F)

    #TODO:
    # - Override true_minimum_distance?

class SymmetricReversibleCode(CyclicCode):
    r"""
    Class for Symmetric Reversible polynomial style cyclic codes.

    INPUT:

    - ``F`` - Field
    - ``n`` - code length
    - ``g=None`` - generator polynomial
    - ``D=None`` - Defining set of the code (might be partially defined)

    Either g or D must be supplied at construction. If D is supplied, the
    actual defining set used is the smallest set possible such that the
    resulting generator matrix fits in F. This is done by taking the cyclotomic
    coset closure of D.
    """
    
    def __init__(self, F, n, g=None, D=None):
        # Check that there exists an m such that n | q^m + 1
        # This is the only requirement for a cyclic code to be symmetric
        # reversible
        m = find_modulo_power(n, 1, F.cardinality())
        if m:
            super(SymmetricReversibleCode, self).__init__(F, n, g=g, D=D)
        else:
            raise Exception("No m such that n | q^m + 1")
