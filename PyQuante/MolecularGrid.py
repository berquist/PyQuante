#!/usr/bin/env python
"""\
 MolecularGrid.py Construct molecular grids from patched atomic
 grids. The technique behind this is based upon:
  A.D. Becke, 'A multicenter numerical integration scheme for
   polyatomic molecules.' J. Chem. Phys 88(4) 1988.

 The atomic grids are constructed from atomic grids that use
 Lebedev grids for the angular part, and Legendre grids for
 the radial parts.

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
from math import sqrt
from AtomicGrid import AtomicGrid, Bragg
from NumWrap import array,concatenate,reshape,zeros
from PyQuante.cints import dist2

class MolecularGrid:
    "Class to hold grid information from patched atomic grids"
    def __init__(self, atoms, nrad=32, fineness=1,**opts):
        self.do_grad_dens = opts.get('do_grad_dens',False)
        self.atoms = atoms
        self.nrad = nrad 
        self.fineness = fineness
        self.make_atom_grids(**opts)
        self.patch_atoms(**opts)
        self._length = None
        return

    def __len__(self):
        self._length = 0
        for agr in self.atomgrids:
            self._length += len(agr)
        return self._length

    def make_atom_grids(self,**opts):
        self.atomgrids = []
        opts['nrad'] = self.nrad
        opts['fineness'] = self.fineness
        for atom in self.atoms:
            atom.grid = AtomicGrid(atom, **opts)
            self.atomgrids.append(atom.grid)
        return

    def patch_atoms(self,**opts):
        do_becke_hetero = opts.get('do_becke_hetero',False)
        nat = len(self.atoms)
        for iat in range(nat):
            ati = self.atoms[iat]
            npts = len(self.atomgrids[iat])
            for i in range(npts):
                point = self.atomgrids[iat].points[i]
                xp,yp,zp,wp = point.xyzw()
                rip2 = dist2(ati.pos(),(xp,yp,zp))
                rip = sqrt(rip2)
                sprod = 1
                for jat in range(nat):
                    if jat == iat: continue
                    atj = self.atoms[jat]
                    rjp2 = dist2(atj.pos(),(xp,yp,zp))
                    rjp = sqrt(rjp2)
                    rij2 = dist2(ati.pos(),atj.pos())
                    rij = sqrt(rij2)
                    mu = (rip-rjp)/rij
                    # Modify mu based on Becke hetero formulas (App A)
                    if do_becke_hetero and ati.atno != atj.atno:
                        chi = Bragg[ati.atno]/Bragg[atj.atno]
                        u = (chi-1.)/(chi+1.)
                        a = u/(u*u-1)
                        a = min(a,0.5)
                        a = max(a,-0.5)
                        mu += a*(1-mu*mu)
                    sprod *= sbecke(mu)
                    #if rjp2 < rip2: point.flag = True
                point._w *= sprod
        return
    

    def points(self):
        "Dynamically form an array of all grid points"
        p = []
        for agr in self.atomgrids: p.extend(agr.points)
        return p

    def set_bf_amps(self,bfs,**opts):
        "Set the basis func amplitude at each grid point"
        for agr in self.atomgrids: agr.set_bf_amps(bfs,**opts)
        return

    def setdens(self,D,**opts):
        "Set the density at each grid point"
        for agr in self.atomgrids: agr.setdens(D,**opts)
        return

    def weights(self):
        "Return a vector of weights of each point in the grid"
        weights = array((),'d')
        for agr in self.atomgrids:
            aw = agr.weights()
            weights = concatenate((weights,aw))
        return weights
    
    def dens(self):
        "Return the density for each point in the grid"
        ds = array((),'d')
        for agr in self.atomgrids:
            ad = agr.dens()
            ds = concatenate((ds,ad))
        return ds

    def gamma(self):
        "Return the density gradient gamma for each point in the grid"
        if not self.do_grad_dens: return None
        gs = array((),'d')
        for agr in self.atomgrids:
            ag = agr.gamma()
            gs = concatenate((gs,ag))
        return gs

    def grad(self):
        pts = self.points()
        npts = len(pts)
        gr = zeros((npts,3),'d')
        for i in range(npts):
            gr[i,:] = pts[i].grad()
        return gr        

    def gradbfab(self,ibf,jbf):
        "Computes grad(ibf*jbf) over the grid"
        assert self.do_grad_dens
        pts = self.points()
        npts = len(pts)
        grab = zeros((npts,3),'d')
        for i in range(npts):
            pt = pts[i]
            grab[i,:] = pt.bfs[ibf]*pt.bfgrads[jbf,:] +\
                        pt.bfs[jbf]*pt.bfgrads[ibf,:]
        return grab

    def bfgrad(self,ibf):
        pts = self.points()
        npts = len(pts)
        gra = zeros((npts,3),'d')
        for i in range(npts):
            gra[i,:] = pts[i].bfgrads[ibf,:]
        return gra        
    
    def bfs(self,i):
        "Return a basis function over the entire grid"
        bfs = array((),'d')
        for agr in self.atomgrids:
            abfs = agr.bfs(i)
            bfs = concatenate((bfs,abfs))
        return bfs

    def nbf(self):
        return self.atomgrids[0].nbf()

    def npts(self):
        npts = 0
        for agr in self.atomgrids: npts += agr.npts()
        return npts

    def allbfs(self):
        "Construct a matrix with bfs in columns over the entire grid, "
        " so that R[0] is the first basis function, R[1] is the second..."
        bfs = array((),'d')
        for agr in self.atomgrids:
            abfs = agr.allbfs()
            bfs = concatenate((bfs,abfs))
        # Now the bfs array is a concatenation of all of the bfs
        npts = self.npts()
        nbf,nrem = divmod(len(bfs),npts)
        if nrem != 0: raise "Remainder in divmod allbfs"
        nbf2 = self.nbf()
        if nbf != nbf2: raise "Wrong # bfns %d %d" % (nbf,nbf2)
        bfs = reshape(bfs,(npts,nbf))
        return bfs

# These are the functions for the becke projection operator
def fbecke(x,n=3):
    for i in range(n): x = pbecke(x)
    return x
def pbecke(x): return 1.5*x-0.5*pow(x,3)
def sbecke(x,n=3): return 0.5*(1-fbecke(x,n))

if __name__ == '__main__':
    # Test the becke projection grids
    from PyQuante.Molecule import Molecule
    h2 = Molecule('h2',
                  atomlist = [(1,(0.,0.,0.7)),(1,(0.,0.,-0.7))],
                  units = 'Bohr')
    grid = MolecularGrid(h2,do_becke=True)
    
