#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Generally usable polynomial functions"""
from codinglib.util import *
from codinglib.ea import *

def decompose(f, (g1, g2)):
    """Find h1, h2 in F[x] such that f = h1*g1 + h2*g2 where F is the base
    field and f, g1, g2 polynomials in x over this field. The h1 and h2 are
    sought so that the maximal of their degrees is minimised"""
    #TODO: Currently just minimises the degree of h1
    eaRes = ea(g1, g2)
    gcd = eaRes[0][0]
    if not gcd.is_one():
        if not (f % gcd).is_one():
            raise Exception("Target polynomial not describable""")
        return decompose(f/gcd, g1/gcd, g2/gcd)


def weighted_degree(self, *weights):
    """
    Return the weighted degree of ``self``, which is the maximum weighted
    degree of all monomials in ``self``; the weighted degree of a monomial
    is the sum of all powers of the variables in the monomial, each power
    multiplied with its respective weight in ``weights''.

    TODO: This implementation is slow. It should be reimplemented to use
    the Singular interface.

    INPUT::

    - ``weights`` - Either individual numbers, a tuple or a dictionary,
      specifying the weights of each variable. If it is a dictionary, it
      maps each variable of ``self'' to its weight. If it is individual
      numbers or a tuple, it specified the weights in the order of the
      generators as given by ``self.parent().gens()'':

    EXAMPLES::

        sage: R.<x,y,z> = GF(7)[]
        sage: p = x^3 + y + x*z^2
        sage: p.weighted_degree({z:0, x:1, y:2})
        3
        sage: p.weighted_degree(1, 2, 0)
        3
        sage: p.weighted_degree((1, 4, 2))
        5
        sage: p.weighted_degree((1, 4, 1))
        4
        sage: q = R.random_element(100, 20) #random
        sage: q.weighted_degree(1,1,1) == q.total_degree() 
        True
    """
    # Make the weights into a vector so we can use dot_product
    if len(weights) ==  1:
        # First unwrap it if it is given as one element argument
        weights = weights[0]

    if isinstance(weights, tuple): 
        weights = vector(weights)
    elif isinstance(weights, dict):
        weights_list = [ weights[g] for g in self.parent().gens() ]
        weights = vector(weights_list)
    else:
        raise Exception('weights must be a tuple or a dictionary')

    # Go through each monomial, calculating the weight by dot product
    deg = max([ weights.dot_product(vector(m.degrees()))
                    for m in self.monomials() ])
    return deg
    

def companion_polynomial(p):
    """
    Return the companion polynomial to p, p' = x^deg(p)*p(1/x)

    INPUT::
        p : A univariate polynomial    

    OUTPUT::
        The companion polynomial to p in the same variable.
    """
    # Coerce p to a multivariate ring in one variable to get to the homogenize
    # method there
    x = p.parent().gen()
    from sage.rings.polynomial.multi_polynomial_libsingular \
        import MPolynomialRing_libsingular
    print "Parent = %s , Ring = %s" % (p.parent(), p.base_ring())
    MPF = MPolynomialRing_libsingular(p.base_ring(), 1, 'z')
    z = MPF.gen()
    q = MPF(p(z)).homogenize()
    nvar = q.nvariables()
    if nvar == 2:
        return q(1,x)
    elif nvar == 1:
        # Happens if p -- and therefore q -- are monomials
        return q(1)
    else:
        # This happens only if p or q is a constant, so its companion is itself
        return p 
