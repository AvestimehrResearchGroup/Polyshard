#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Class for representing bipartite graph codes with corresponding utility
functions"""
from codinglib.util import *
from codinglib.code import *
import codinglib.codeTesting

class GraphCode(BlockCodeAbstract):
    r"""
    Class for graph code constructions: Given a bipartite graph G, where each
    left-vertex has degree d1 and each right-vertex has degree d2, a left
    component code C1 of length d1 and a right component code C2 of length d2,
    both codes over the same field F, consider placing field elements on each
    edge such that the edges incident to some left-vertex is a codeword of C1,
    while the edges incident to some right-vertex is a codeword of C2.
    
    INPUT:

    - ``G`` - The bipartite graph, given as an incidence matrix n1xn2, where
              the (i,j) is 1 iff there is an edge between left-vertex i and
              rightvertex j. This naturally induces the ordering of the
              vertices, and we then canonically use the ordering of edges such
              that edge (l1, r1) < (l2, r2) iff l1 < l2 or l1 = l2 and r1 < r2.
    - ``C1`` - The left component code.
    - ``C2`` - The right component code.

    FIELDS:

    - ``F`` - Field
    - ``n`` - Code length
    - ``G`` - The bipartite graph
    - ``C1`` - The left component code.
    - ``C2`` - The right component code.
    - ``edgemap`` - Given a (l1, r1) which is an edge in the graph, returns the
                    number this edge has been given in the canonical ordering.
    """
    def __init__(self, G, C1, C2):
        self.G = G
        self.C1 = C1
        self.C2 = C2

        assert(G.base_ring() == ZZ)
        for r in G.rows():
            assert(all(ri == 0 or ri == 1 for ri in r))

        n1, n2 = C1.n, C2.n
        self.n = n1*G.nrows()
        assert( self.n == n2*G.ncols() )

        self.F = C1.F
        assert( self.F==C2.F )

        def checkdegs(Gi, d):
            degs = Gi * vector([1]*Gi.ncols())
            if not all(di == d for di in degs):
                raise Exception("Incidence degree should be the same for all"
                    + " vertices of the same side")
        checkdegs(G, n1)
        checkdegs(G.transpose(), n2)

        edgemap = dict()
        i = 0
        for r in range(0,G.nrows()):
            for c in range(0,G.ncols()):
                if G[r, c].is_one():
                    edgemap[(r, c)] = i
                    i = i+1
        self.edgemap = edgemap

    def __str__(self):
        return ("[%s,?,?] bipartite craph code with left component code\n\t%s "\
              +"\nand right component code\n\t%s") % (self.n,self.C1,self.C2)

    @cached_method
    def parity_check_matrix(self):
        """Return the parity check matrix for the graph code."""
        H1 = self.C1.parity_check_matrix()
        H2 = self.C2.parity_check_matrix()
        p1, p2 = H1.nrows(), H2.nrows()
        G = self.G
        n1 = G.nrows()
        n2 = G.ncols()
        edgemap = self.edgemap

        part1rows = p1*n1
        bigH = matrix(self.F,part1rows + p2*n2,self.n)

        #Add row parity checks
        for r in range(0, n1):
            myH = matrix(self.F, p1, self.n)
            hc = 0
            for c in range(0, n2):
                if G[r,c].is_one():
                    myH.set_column(edgemap[r,c], H1.column(hc))
                    hc = hc+1
            bigH.set_block(r*p1,0,myH)

        #Add column parity checks
        for c in range(0, n2):
            myH = matrix(self.F, p2, self.n)
            hc = 0
            for r in range(0, n1):
                if G[r,c].is_one():
                    myH.set_column(edgemap[r,c], H2.column(hc))
                    hc = hc+1
            bigH.set_block(part1rows+c*p2,0,myH)

        return bigH
            
    def true_dimension(self):
        """Return the true dimension of the code by calculating the parity
        check matrix and deducting its rank from the code length"""
        return self.n - self.parity_check_matrix().rank()
        





