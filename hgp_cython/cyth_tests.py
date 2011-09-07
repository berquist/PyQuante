# Can I speed up the vrr routines using cython?

from itertools import combinations_with_replacement,product,izip
from PyQuante import hgp,chgp
import hgp_cython
import time

def approx(a,b,delta=1e-7): return abs(a-b)<delta
def same(l,delta=1e-7):
    return False not in [approx(l[i],l[0],delta) for i in range(1,len(l))]

def test_match(*args):
    n = len(args[0])
    nargs= len(args)
    print "comparing %d methods, each with %d integrals" % (nargs,n)
    nerr = 0
    for inti in izip(*args):
        if not same(inti):
            nerr += 1
    print "total of %d errors found" % nerr
    return

def allints(method):
    xa,ya,za = 0.1,0.2,0.3
    xb,yb,zb = 0.3,0.4,0.5
    xc,yc,zc = 0.4,0.5,0.6
    xd,yd,zd = 0.5,0.6,0.7
    norma=normb=normc=normd=1.0
    alphaa,alphab,alphac,alphad = 1.,2.,3.,4.
    ints = []
    M = 0
    for la,ma,na in product(xrange(3),repeat=3):
        for lc,mc,nc in product(xrange(3),repeat=3):
            ints.append(method.vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                   (xd,yd,zd),normd,alphad,M))
    return ints

t0 = time.time()
ints_hgp = allints(hgp)
t1 = time.time()
ints_chgp = allints(chgp)
t2 = time.time()
ints_cython = allints(hgp_cython)
t3 = time.time()
test_match(ints_hgp,ints_chgp,ints_cython)
print "Timing: hgp %.4f chgp %.4f cython %.4f" % (t1-t0,t2-t1,t3-t2)
