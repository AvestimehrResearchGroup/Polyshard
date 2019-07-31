# Test that CountingFields behave exactly like a Field, and that we never
# inadvertently escape from the counting field

F = GF(7)
CF = CountingField(F)
a = CF.gen()
assert(isinstance(a, CountingFieldElement))
assert(isinstance(CF.random_element(), CountingFieldElement))

#
CF2 = CountingField(F)
assert(CF==CF2)
#

#PF.<x> = CF[]
PF_d  = PolynomialRing(CF, 'x')
# assert (PF.base_ring() == CF and PF_d.base_ring() == CF)

#

a==a
a+a
a*a

#

a-a
a/a

#

assert( isinstance(CF(0), CountingFieldElement) )
assert( isinstance(CF(1), CountingFieldElement) )

#

print CF.zero_element() == 0
print CF(0) == 0
print a-a==0
print a+0==a
print a*0==0
print a*1==a
print a*(-1)==-a
print -a==0-a
print 0/a==0

### TODO: Test that coercion from pure field fails

### Test polynomials over CF

PF.<x> = CF[]
p = PF(a + x)
p = PF(x + a)
p = PF.random_element()

print p + CF(0)
print p + 0
print p + 1
print p(a)
print 
print p*(2*p)
print 
print p.roots(multiplicities=False)

### Test matrices over CF

M = matrix(CF, 3, 3, lambda i,j: CF.random_element() )
M*M
M+3*M

### Test matrices over CF[]

M = matrix(PF, 3, 3, lambda i,j: PF.random_element(degree=4) )
M*M
(x+a)*M+a*M


### Plot field operations used for polynomial multiplication or other ops
# These tests demonstrate that Sage internally computes multiplication of
# polynomials FAST but division SLOW.
import time
before = time.time()
def op(p):
    #p*p
    p2 = PF.random_element(degree=2*p.degree())
    p.base_ring().reset_counters()
    p2 % p
CF = CountingField(GF(7))
PF.<x> = CF[]
def get_poly_of_deg(deg):
    p = PF.random_element(degree=deg)
    if p.degree() < deg:
        return p + x^deg
    else:
        return p
def count_multiplication(deg):
    CF.reset_counters()
    p = get_poly_of_deg(deg)
    op(p)
    return CF.multiplications()

var('d')
maxdeg = 100
g = int_plot(count_multiplication, ('d', 1, maxdeg))
g_bs = g + plot(.5*d^2, (d,0,maxdeg), color='red')
# g_bs = g + plot(4.5*d*log(d), (d,0,maxdeg), color='red')
# g_bs += plot(7*d*log(d), (d,0,maxdeg), color='red')
print "Took %f seconds" % (time.time() - before)




### Plot field operations used for mass-evaluating polynomials over various fields
import time
before = time.time()
degs = [100]
qmax = 500
qmin = 50
def get_poly_of_deg(deg, PF):
    (x,) = PF.gens()
    p = PF.random_element(degree=deg)
    if p.degree() < deg:
        return p + x^deg
    else:
        return p
counts = dict()
for N in range(qmin, qmax+1):
    try:
        F = CountingField(GF(N))
        PF.<x> = F[]
        M = poly_evaluate_build_M(F.list(), PF)
        for deg in degs:
            p = get_poly_of_deg(deg, PF)
            p = p.monic()
            F.reset_counters()
            evs_n = poly_evaluate(p, F.list(), algorithm='naive')
            mlts_n = F.multiplications()
            F.reset_counters()
            evs_gg = poly_evaluate(p, F.list(), algorithm='GG', M=M)
            mlts_gg = F.multiplications()
            assert evs_n == evs_gg
            counts[(deg,N)] = (mlts_n, mlts_gg)
    except ValueError: #no field of this size
        pass
# Plot results
# Note: the zero values on the plot correspond to values we didn't test
var('N')
g = Graphics()
for i in range(len(degs)):
    deg = degs[i]
    g += int_plot(lambda N: counts[(deg,N)][0] if (deg,N) in counts else 0, (N, qmin, qmax), color=util.some_colors[i], legend_label="deg = %s" % degs[i])
    g += int_plot(lambda N: counts[(deg,N)][1] if (deg,N) in counts else 0, (N, qmin, qmax), color=util.some_colors[i], legend_label="deg = %s" % degs[i], marker='D')
    g += plot(N*(degs[i]-1), (N,qmin, qmax), color=util.some_colors[i])
print "Took %f seconds" % (time.time() - before)



### Compare PolynomialRing.lagrange_polynomial with manual method using trivial caching
F = CountingField(GF(53, 'alpha'))
PF.<x> = F[]
word = [ F.random_element() for i in range(F.cardinality()) ]
F.reset_counters()
Fl = F.list()
R = PF.lagrange_polynomial( zip(Fl, word) )
lagrange = F.reset_counters()[1]
Hs = [ prod(x - b for b in F.list() if b != a) for a in F.list() ]
etas = [ Hs[i](Fl[i]) for i in range(F.cardinality()) ]
Hns = [ Hs[i]/etas[i] for i in range(F.cardinality()) ]
F.reset_counters()
R2 = sum( word[i]*Hns[i] for i in range(F.cardinality()) )
manual = F.multiplications()
print "Lagrange: %s, Manual: %s" % (lagrange, manual)





### Perform full Gao decoding
F = CountingField(GF(53, 'alpha'))
PF.<x> = F[]
n, k = 15, 7
C = GRS(F, n, k, F.list()[1:n+1])
tau = (n-k)//2

G = prod(x-alpha for alpha in C.alphas) # Offline, don't count

PF.<x> = F[]
def get_decoding_instance():
    c = C.random_codeword()
    e = random_error(n,F,tau)
    r = c + e
    return c, e, r


c, e, r = get_decoding_instance()
F.reset_counters()
R = PF.lagrange_polynomial(zip(C.alphas, r))
Mw = matrix(PF, 2, 2, [[x^k, R],[0,G]])
module_weak_popov(Mw)

i = 0 if LP(Mw.row(0))==0 else 1
f = PF( Mw[i, 1] // Mw[i, 0].shift(-k) )

print "Gao Multiplications: ", F.multiplications()
print "Gao Additions: ", F.additions()
assert( C.encode(poly2list(f, C.k)) == c )



### Perform Gao decoding with re-encoding
# Some offline computations, don't count
G_r = prod(x-alpha for alpha in C.alphas[k:])
g = prod(x-alpha for alpha in C.alphas[:k])
hs_r = [ G_r//(x-a) for a in C.alphas[k:] ]
Hs_r = [ hs_r[i-k]/(hs_r[i-k](C.alphas[i]))/g(C.alphas[i]) for i in range(k,n) ]

c, e, r = get_decoding_instance()
F.reset_counters()

# Reencode
f_r = PF.lagrange_polynomial(zip(C.alphas[:k], r[:k]))
r_r = vector([ F.zero() ]*k + [ r[i] - f_r(C.alphas[i]) for i in range(k,n) ])

R_r = sum( r_r[i] * Hs_r[i-k] for i in range(k,n) )
Mw = matrix(PF, 2, 2, [[x*R_r, 1],[x*G_r, 0]])
module_weak_popov(Mw)

j = 0 if LP(Mw.row(0))==1 else 1
Lambda = Mw[j, 1]

print "Gao Re-encoded Multiplications: ", F.multiplications()
print "Gao Re-encoded Additions: ", F.additions()
assert( set([ i for i in range(n) if Lambda(C.alphas[i]).is_zero() ]) == support(e) )



### Perform Syndrome decoding using Berlekamp--Massey
c, e, r = get_decoding_instance()
C.parity_colmults() # cache-calculate parity column multipliers
F.reset_counters()

S = C.syndrome(r, direct=False)
(Lambda_seq, L, _) = bma(S)

print "BMA Multiplications (no root finding): ", F.multiplications()
print "BMA Addition (no root finding): ", F.additions()
Lambda_seq.reverse()
Lambda = PF(Lambda_seq)

assert( set([ i for i in range(n) if Lambda(C.alphas[i]).is_zero() ]) == support(e) )


### Perform Guruswami-Sudan decoding
c, e, r = get_decoding_instance()
D = GRSDecoderGuruswamiSudan(C, tau=5)
F.reset_counters()

out = D.decode(r)

print "G-S(s = %s) Multiplications: %s" % (D.s, F.multiplications())
print "G-S(s = %s) Additions: %s" % (D.s, F.additions())
assert c in out


### Perform Power decoding
c, e, r = get_decoding_instance()
D = GRSDecoderPower(C, gs_params(n, k, 5))
F.reset_counters()

cout = D.decode(r)

print "PowerGao(s = %s) Multiplications: %s" % (D.s, F.multiplications())
print "PowerGao(s = %s) Additions: %s" % (D.s, F.additions())
assert c == cout


### Perform PowerSyndrome decoding
c, e, r = get_decoding_instance()
D = GRSDecoderPowerSyndrome(C, gs_params(n, k, 5))
F.reset_counters()

cout = D.decode(r)

print "PowerSyndrome(s = %s) Multiplications: %s" % (D.s, F.multiplications())
print "PowerSyndrome(s = %s) Additions: %s" % (D.s, F.additions())
assert c == cout

