# Test the two electron integral methods to insure that they all
# produce the same results for a variety of cases
from itertools import combinations_with_replacement
from PyQuante import pyints,hgp,rys,cints,chgp,crys
from PyQuante.PGBF import PGBF

def approx(a,b,delta=1e-7): return abs(a-b)<delta
def same(args):
    return False not in [approx(args[i],args[0]) for i in range(len(args))]

def test(a,b,c,d,sa,sb,sc,sd,ints):
    data = []
    for lib in int_libs:
        datum = lib.coulomb_repulsion(a.origin,a.norm,a.powers,a.exp,
                                      b.origin,b.norm,b.powers,b.exp,
                                      c.origin,c.norm,c.powers,c.exp,
                                      d.origin,d.norm,d.powers,d.exp)
        data.append(datum)
    if not same(data):
        print sa,sb,sc,sd,data
    return 

def test_bfs(bfs,syms,ints):
    assert len(bfs) == len(syms)
    for i,j,k,l in combinations_with_replacement(xrange(len(bfs)),4):
        test(bfs[i],bfs[j],bfs[k],bfs[l],syms[i],syms[j],syms[k],syms[l],ints)
    return

# Basis functions over which to test. The random-ish numbers here just insure
# that I'm not limiting myself to certain symmetry constraints.
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
xxx = PGBF(0.35,(0.15,0.15,0.15),(3,0,0))
yyy = PGBF(0.25,(0.2,0.2,0.2),(0,3,0))
zzz = PGBF(0.15,(0.1,0.1,0.1),(0,0,3))
xyz = PGBF(0.5,(0.2,0.2,0.2),(1,1,1))

syms = ['s1','s2','s3','px','py','pz','xx','yy','zz','xy','yz','xz']
fsyms = ['xxx','yyy','zzz','xyz']

bfs = [s1,s2,s3,px,py,pz,xx,yy,zz,xy,yz,xz]
fbfs = [xxx,yyy,zzz,xyz]

int_libs = [pyints,hgp,cints,chgp]
rys_libs = [rys,crys]

print "testing all libraries for non-f-functions"
test_bfs(bfs,syms,int_libs+rys_libs)
print "just testing non-rys libraries for basis set that contains f-functions"
test_bfs(bfs+fbfs,syms+fsyms,int_libs)

#print "testing permutations of the zzz,xyz case that seems to be the worst"
#test(zzz,zzz,xyz,xyz,'zzz','zzz','xyz','xyz',int_libs)
#test(xyz,xyz,zzz,zzz,'xyz','xyz','zzz','zzz',int_libs)
#test(zzz,xyz,zzz,xyz,'zzz','xyz','zzz','xyz',int_libs)
#test(xyz,zzz,xyz,zzz,'xyz','zzz','xyz','zzz',int_libs)
