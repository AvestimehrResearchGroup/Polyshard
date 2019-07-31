### Test GRS code implementation

### Test equality
F = GF(29)
n, k = 25, 19
C1 = GRS(F, n, k, F.list()[:n], [F(i) for i in range(n)])
F2 = GF(29)
C2 = GRS(F2, n, k, F2.list()[:n], [F2(i) for i in range(n)])

assert C1==C2, "eq 1"

C3 = GRS(F, n, k-1, F.list()[:n], [F2(i) for i in range(n)])
assert C1!=C3 , "neq 3"
C4 = GRS(F, n-1, k, F.list()[:n-1], [F2(i) for i in range(n-1)])
assert C1!=C4 , "neq 4"
C5 = GRS(F, n, k, F.list()[:n])
assert C1!=C5 , "neq 5"
C6 = GRS(F, n, k, F.list()[1:n+1], [F2(i) for i in range(n)])
assert C1!=C6 , "neq 6"
