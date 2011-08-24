# Test the two electron integral methods
from itertools import combinations_with_replacement
from PyQuante import pyints,hgp,rys,cints,chgp,crys
from PyQuante.PGBF import PGBF

def approx(a,b,delta=1e-7): return abs(a-b)<delta

s1 = PGBF(1.0,(0,0,0),(0,0,0))
s2 = PGBF(2.0,(.1,.2,.3),(0,0,0))
s3 = PGBF(3.0,(.1,.2,.3),(0,0,0))
px = PGBF(3.0,(.2,.4,.5),(1,0,0))
py = PGBF(4.0,(.1,.5,.9),(0,1,0))
pz = PGBF(5.0,(.3,.2,.1),(0,0,1))
xx = PGBF(1.1,(.1,.2,.3),(2,0,0))
yy = PGBF(1.3,(1.1,.2,.3),(0,2,0))
zz = PGBF(1.5,(.1,2.2,.3),(0,0,2))
xy = PGBF(1.7,(.1,1.2,4.3),(1,1,0))
yz = PGBF(1.2,(.1,1.2,4.3),(0,1,1))
xz = PGBF(1.4,(3.1,5.2,.3),(1,0,1))
syms = ['s1','s2','s3','px','py','pz','xx','yy','zz','xy','yz','xz']
bfs = [s1,s2,s3,px,py,pz,xx,yy,zz,xy,yz,xz]
nbf = len(bfs)
for i,j,k,l in combinations_with_replacement(xrange(nbf),4):
    a,b,c,d = bfs[i],bfs[j],bfs[k],bfs[l]
    data = []
    for lib in [pyints,hgp,rys,cints,chgp,crys]:
        datum = lib.coulomb_repulsion(a.origin,a.norm,a.powers,a.exp,
                                   b.origin,b.norm,b.powers,b.exp,
                                   c.origin,c.norm,c.powers,c.exp,
                                   d.origin,d.norm,d.powers,d.exp)
        data.append(datum)
    i1,i2,i3,i4,i5,i6 = data
    if not (approx(i1,i2) and approx(i1,i3) and approx(i1,i4) and
            approx(i1,i5) and approx(i1,i6)):
        print syms[i],syms[j],syms[k],syms[l],i1,i2,i3,i4,i5,i6

