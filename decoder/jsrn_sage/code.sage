#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Abstract and concrete class for representing linear block codes"""
from codinglib.util import *

#TODO: Make central __str__, __repr__ and _latex_ by having a field code_class_name

class BlockCodeAbstract(object):
    # Abstract class for general linear block codes. This is the main
    # class for sub-classing when implementing new families of codes.
    #
    # A certain amount of processing is expected upon construction of a code
    # class, but the idea is (as in many instances in Sage), that construction
    # itself of a Code should be cheap, and as much computation as possible
    # deferred to when one actually calls its functions. For instance, avoid
    # computing the generator matrix of an indirectly defined (e.g. algebraic)
    # code in the constructor. Instead, use @cached_method freely to avoid
    # multiple computations of the same object.
    #
    # Don't construct codes using this class: perhaps you are looking
    # for BlockCode
    #
    # Fields
    # F : field of elements
    # n : code length
    # k : designed code dimension
    # d : designed minimum distance
    #
    # The following methods are possible to overwrite
    # - generator_matrix *
    # - parity_check_matrix *
    # - encode
    # - syndrome
    # - iscodeword
    # - true_dimension *
    # - true_minimum_distance *
    # - random_codeword
    #
    # For those methods with stars, it is important that they are defined with
    # @cached_method to avoid multiple calculations of the same object.
    #
    # The following decribes which methods standard implementations fall back
    # to if they are not present
    # - generator_matrix: parity_check_matrix
    # - parity_check_matrix: generator_matrix
    # - encode: generator_matrix
    # - iscodeword: syndrome
    # - random_codeword: generator_matrix


    def __init__(self):
        pass

    def __str__(self):
        raise NotImplementedError

    def __contains__(self, x):
        return self.iscodeword(x)

    def __eq__(self, other):
        """Very superficial, but relatively cheap, check for total equivalence of two codes: their generator matrices are exactly the same."""
        return self.generator_matrix() == other.generator_matrix()

    @cached_method
    def generator_matrix(self):
        """Return a parity check matrix for the code."""
        if hasattr(self,"defaulting_parity_check"):
            raise NotImplementedError
        else:
            return self.parity_check_matrix().right_kernel().basis_matrix()

    @cached_method
    def parity_check_matrix(self):
        """Return a generator matrix for the code."""
        self.defaulting_parity_check = True
        Sp = self.generator_matrix().right_kernel()
        return Sp.basis_matrix()

    def encode(self, info):
        """Encode the information word"""
        G = self.generator_matrix()
        return vector(info)*G

    @cached_method
    def _an_information_set(self):
        """Return an information set for G. Repeated calls to this method will
        return the same information set."""
        #Randomised, optimistic algorithm: draw k random columns and ask if they
        #have full rank.
        k = self.true_dimension()
        pos = list(range(self.n))
        G = self.generator_matrix()
        while True:
            shuffle(pos)
            info_set = pos[:k]
            info_set.sort()
            Gt = G.matrix_from_columns(info_set)
            if rank(Gt) == k:
                return info_set

    @cached_method
    def _unencoder_matrix(self):
        """Find an information set for G, and return the inverse of those
        columns of G: this allows o unencode codewords"""
        Gt = self.generator_matrix().matrix_from_columns(self._an_information_set())
        return Gt.inverse()
        
    def unencode(self, c):
        """Return the information corresponding to a codeword. This function
        satisfies `unencode(encode(info)) == info`.
        Since the information word is calculated from an information set of the
        code, then if `c` is not a codeword, then the returned vector might not
        be closer than `n-k` to `c`"""
        U = self._unencoder_matrix()
        info_set = self._an_information_set()
        cc = vector( c[i] for i in info_set )
        return cc * U

        #TODO: Rename to is_codeword
    def iscodeword(self, r):
        """Returns whether r is a codeword"""
        # Check that the word is in the field
        if not all(ri in self.F for ri in r):
            return false
        # Check that the syndrome is zero
        syn = self.syndrome(r)
        if hasattr(syn, "is_zero"):
            return syn.is_zero()
        elif isinstance(syn, list):
            return all(s.is_zero() for s in syn)
        else: 
            raise Exception("Don't know how to check if syndrome is zero.")

    def syndrome(self, r):
        """Returns the syndrome of the received word.
        INPUT:
        - ``r`` -- The word to be evaluated, in the preferred form.

        OUTPUT:
        The syndrome in a preferred form. If the standard implementation of
        iscodeword is to be used, it should be a an object with the is_zero
        function, or a list of such"""
        return self.parity_check_matrix() * r

    @cached_method
    def true_dimension(self):
        """Return the true dimension of the code by calculating the generator
        matrix and calculating its rank"""
        return self.generator_matrix().rank()

    def random_codeword(self):
        """Return a random codeword from the code. Requires the calculation of
        the generator matrix"""
        k = self.true_dimension();
        return random_vector(self.F, k) * self.generator_matrix();

    @cached_method
    def true_minimum_distance(self):
        """Compute the true minimum distance of this code by analysing the
        generator matrix. Computation time is exponential in code dimension."""
        # Calculate by calling the Sage integrated functions
        Csage = LinearCode(self.generator_matrix())
        return Csage.minimum_distance()

    def list(self):
        """Tabulate and return all codewords of the code.
        NOTE: This list will have len (F.cardinality())^k"""
        G = self.generator_matrix()
        RS = G.row_space()
        return RS.list()
        
    def __iter__(self):
        for c in self.list():
            yield c


class BlockCode(BlockCodeAbstract):
    r"""
    Class for general linear block codes construced using a generator matrix or parity check matrix
    Very thin layer on top of Code

    INPUT:

    - ``M`` - A matrix which will be either the generator or parity check
    - ``generator_matrix=True`` - Indicates whether ``M`` is a generator or parity check matrix

    FIELDS:
    
    - F - Field of elements
    - n - Code length
    - k - Designed code dimension (determined directly by the number of rows in M)
    """

    def __init__(self, M, generator_matrix=True):
        """Construct a code using the given matrix. If
        `generator_matrix==True` then `M` is interpreted as the
        generator matrix of the code, otherwise the parity check
        matrix"""
        if generator_matrix:
            self.G = M
            self.n = M.ncols()
            self.k = M.nrows()
        else:
            self.H = M
            self.n = M.ncols()
            self.k = self.n - M.nrows()
        self.F = M.base_ring()
        self.d = None

    def __str__(self):
        d = self.d if self.d else "??"
        return "[%s,>=%s,%s] block code over %s" % (self.n,self.k,d,self.F)

    def __repr__(self):
        return self.__str__()

    def _latex_(self):
        d = self.d if self.d else "??"
        return r"[%s,\ge %s,%s]{\rm\ block \ code\ over\ } %s" % (self.n,self.k,d,latex(self.F))


    @cached_method
    def generator_matrix(self):
        """Return a generator matrix for the code."""
        try:
            return self.G
        except:
            return self.parity_check_matrix().right_kernel().basis_matrix()

    @cached_method
    def parity_check_matrix(self):
        """Return a parity check matrix for the code."""
        try:
            return self.H
        except:
            Sp = self.generator_matrix().right_kernel()
            return Sp.basis_matrix()

    
class DecodingFailedError(Exception):
    pass

class DecoderBlockCodeAbstract:
    # Abstract class for decoding algorithms for linear block codes. This is the
    # main class for sub-classing when implementing new decoding algorithms for
    # codes.
    #
    # There are a number of advantages for making decoding algorithms into
    # classes:
    # 1) It allows for preprocessing when a given decoder is applied repeatedly.
    # 2) It allows for giving additional decoder arguments at construction time,
    # letting the decode function itself having uniform calling convention.
    # 3) It allows to add properties to decoding algorithms which can be queried
    # (e.g. is_list_decoder).
    #
    # Fields
    # C : code for which the decoder works
    #
    # The following methods are possible to overwrite
    # - decode
    # - decode_to_information
    # - is_list_decoder
    #
    # The following decribes which methods standard implementations fall back
    # to if they are not present
    # - decode: decode_to_information
    # - decode_to_information: decode

    def __init__(self):
        pass

    def __str__(self):
        return "Decoding algorithm for " + str(self.C)

    def decode(self, r):
        """Decode the received word r into a codeword from C (or a list of such,
        if is_list_decoder is true). May raise DecodingFailedError"""
        if hasattr(self,"_defaulting_decode"):
            raise NotImplementedError
        else:
            if self.is_list_decoder():
                ls = self.decode_to_information(r)
                if ls is None:
                    raise DecodingFailedError("No close codewords found")
                else:
                    return [ self.C.encode(inf) for inf in ls ]
            else:
                return self.C.encode(self.decode_to_information(r))

    def decode_to_information(self, r):
        """Decode the received word r directly into information, such that
        `self.C.encode(inf)` is a codeword close to `r` (or a list of such,
        if is_list_decoder is true). May raise DecodingFailedError"""
        self._defaulting_decode=True
        if self.is_list_decoder():
            ls = self.decode(r)
            return [ self.C.unencode(c) for c in ls ]
        else:
            return self.C.unencode(self.decode(r))
