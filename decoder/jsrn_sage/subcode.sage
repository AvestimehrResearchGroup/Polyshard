#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Direct construction of subfield subcodes by expansion of parity check matrices"""
from codinglib.util import *
from codinglib.code import BlockCodeAbstract

class SubfieldSubcode(BlockCodeAbstract):
    # Abstract class for general linear block codes
    # Don't construct codes using this class: instead sub-class it.
    # Perhaps you are looking for BlockCode
    #
    # Fields
    # Cbig : the super code for which this code is defined
    # F : field of elements of the subfield. Subfield of Cbig.F
    # n : code length
    # d : designed minimum distance, i.e. equals Cbig.d
    #

    def __init__(self, Cbig, F):
        self.Cbig = Cbig
        self.F = F
        if self.F != Cbig.F.prime_subfield():
            raise NotImplementedError("Construction of subfield subcodes which are not over the prime field of the super field")
        H = matrix_subfield(Cbig.parity_check_matrix())
        H = H.rref() # TODO: Inplace would be better
        self._parity_check_matrix = H

    def parity_check_matrix(self):
        return self._parity_check_matrix
