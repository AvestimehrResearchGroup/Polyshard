#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

r"""Chinese remainder codes: poly-alphabetic codes over `ZZ/pZZ` for different
primes `p`, which are the `ZZ`-equivalent of Reed-Solomon codes. Note that this
code implementation does not inherit from BlockCodeAbstract since it is
poly-alphabetic. Also implementation of InterleavedChineseRemainderCode"""

import codinglib.codeTesting

class ChineseRemainderCode:
    r"""Definition of a Chinese Remainder Code. Since it's a poly-alphabetic
    code, it doesn't fit into our structure of Code (so its not a subclass of
    Code), and it doesn't have a 'generator matrix' per se. The codes can be
    classical with $K = \prod(p_i)$ running from $i=0,\ldots,k$ for some $k$, or
    they can be non-classical with any specified $K < N$. In the latter case,
    $k$ will be set to the real number such that $K = (\prod_{i=1}^{k'} P_i )
    P_{k'+1}^{k-k'}$, where $k' = \lfloor k \rfloor$.
    Fields:
        P : The list of moduli
        n : length of P
        N : Product of all p_i in P
        k : The 'dimension' of the code.
        K : The cardinality of the code.
    """
    
    def __init__(self, P, k=None, K=None):
        self.P = copy(P)
        self.n = len(P)
        self.N = prod( p for p in P)
        if k:
            self.k = k
            assert(k < self.n)
            self.K = prod( P[i] for i in range(0,k) )
        elif K:
            self.K = K
            assert(K < self.N)
            # calculate k
            kp=0
            Pp = 1
            while K > Pp:
                Pp *= P[kp]
                kp += 1
            kp -= 1
            Pp = Pp/P[kp]
            if K == Pp:
                self.k = kp
            else:
                self.k = kp + float(log(K/Pp)/log(P[kp]))
        
    def random_codeword(self):
        c = randint(0, self.K-1)
        return vector( c % self.P[i] for i in range(0, self.n) )
        
    def random_error(self, nerrs):
        """Generate a random error with weight nerrs"""
        err_pos = set( codinglib.codeTesting.random_error_pos(self.n, nerrs) )
        return vector( randint(0,self.P[i]-1) if i in err_pos else 0 for i in range(self.n) )
        
    def syndrome(self, r):
        return crt(list(r), self.P)
     
    def information(self, c):
        """Given a codeword, return the information word. Behaviour is unspecified if c is not a codeword"""
        k = crt(list(c), self.P)
        assert( k < self.K or "c was not a codeword")
        return k
    
    def error_locator(self, e):
        """Given an error, return the error locator"""
        return prod( 1 if e[i].is_zero() else self.P[i] for i in range(0, self.n) ) 
               
    def minimum_distance_decoding_bound(self):
        """Return the simple bound for the amount of errors for minimum distance decoder"""
        return floor( (self.n-self.k)*log(self.P[self.k])/( log(self.P[self.k])+log(self.P[-1]) ) )

    def __str__(self):
        return "[%s,%s] CR code with Primes=[%s,...,%s]" % (self.n,self.k, self.P[0], self.P[-1])

    def __repr__(self):
        return self.__str__()

    def _latex_(self):
        return r"[%s,%s]{\rm\ CR\ code\ with\ Primes}$=[%s,...,%s]$" % (self.n,self.k, self.P[0], self.P[-1])



class InterleavedChineseRemainderCode:
    """An interleaved Chinese Remainder Code with same prime set (but possibly different lengths)
    Fields:
        P : The list of moduli
        n : length of P
        N : Product of all p_i in P
        ell : Number of interleaved codes
        ks : List of (non-decreasing) k, each determining cardinality of a
             code. Note that this will be non-integers if the CR codes are not
             classical.
        Ks : List of the cardinalities
        Cs : List of the individual codes
    """
     
    def __init__(self, P, ks=None, Ks=None):
        self.P = copy(P)
        if ks:
            self.ks = copy(ks)
            assert( ks[i] <= ks[i+1] for i in range(len(ks)-1) )
            self.Cs = [ ChineseRemainderCode(P, k=k) for k in ks ]
            self.Ks = [ C.K for C in self.Cs ]
        elif Ks:
            self.Ks = copy(Ks)
            assert( Ks[i] <= Ks[i+1] for i in range(len(Ks)-1) )
            self.Cs = [ ChineseRemainderCode(P, K=K) for K in Ks ]
            self.ks = [ C.k for C in self.Cs]
        else:
            raise Exception("Supply either ks or Ks")
        self.ell = len(self.ks)
        self.n = self.Cs[0].n
        self.N = self.Cs[0].N
    
    def random_codeword(self):
        return matrix(ZZ, [C.random_codeword() for C in self.Cs] )
    
    def syndrome(self, r):
        return vector( self.Cs[i].syndrome( r.row(i) ) for i in range(0,self.ell) )
        
    def information(self, c):
        return vector( self.Cs[i].information( list(c.row(i)) ) for i in range(0, self.ell) )
        
    def random_error(self, nerrs):
        """Generate a random error with weight nerrs"""
        err_cols = set( codinglib.codeTesting.random_error_pos(self.n, nerrs) )
        return matrix(ZZ, self.ell, self.n, lambda i,j: randint(0, self.P[j]-1)  if j in err_cols else 0)
        
    def error_locator(self, e):
        """Given an error, return the error locator"""
        is_err = [ 0 if all( el==0 for el in e.column(i) ) else 1 for i in range(0,self.n) ]
        return prod( 1 if is_err[i]==0 else self.P[i] for i in range(0, self.n) )
        
    def minimum_distance_decoding_bound(self):
        return self.Cs[0].minimum_distance_decoding_bound()


    def __str__(self):
        return "[%s; %s] ICR code with Primes=[%s,...,%s]" % (self.n, ", ".join([str(C.k) for C in self.Cs]), self.P[0], self.P[-1])

    def __repr__(self):
        return self.__str__()

    def _latex_(self):
        return r"[%s; %s]{\rm \ ICR\ code\ with\ Primes}$=[%s,...,%s]$" % (self.n, ", ".join([str(C.k) for C in self.Cs]), self.P[0], self.P[-1])
