from PyQuante import hgp,chgp
from itertools import product

def approx(a,b,delta=1e-7): return abs(a-b)<delta

xa,ya,za = 0.1,0.2,0.3
xb,yb,zb = 0.3,0.4,0.5
xc,yc,zc = 0.4,0.5,0.6
xd,yd,zd = 0.5,0.6,0.7
norma=normb=normc=normd=1.0
alphaa,alphab,alphac,alphad = 1.,2.,3.,4.
MaxL = 5
M = 0
triples = [(6,0,0),(0,6,0),(0,0,6),(5,0,0),(0,5,0),(0,0,5),(4,0,0),(0,4,0),(0,0,4),
           (3,0,0),(0,3,0),(0,0,3),(2,0,0),(0,2,0),(0,0,2),
           (0,0,0),(1,1,1),(2,2,2),(0,0,1),(0,0,2),(0,1,0),(0,2,0)]
nfailed = 0
for la,ma,na in triples:
    #for lc,mc,nc in triples:
    for lc,mc,nc in product(xrange(3),repeat=3):
        for M in [0]:
            py = hgp.vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                         (xb,yb,zb),normb,alphab,
                         (xc,yc,zc),normc,(lc,mc,nc),alphac,
                         (xd,yd,zd),normd,alphad,M)
            c = chgp.vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                         (xb,yb,zb),normb,alphab,
                         (xc,yc,zc),normc,(lc,mc,nc),alphac,
                         (xd,yd,zd),normd,alphad,M)
            if not approx(py,c,1e-12):
                print (la,ma,na),(lc,mc,nc),M,py,c
                nfailed += 1
print "Total of %d cases failed" % nfailed
