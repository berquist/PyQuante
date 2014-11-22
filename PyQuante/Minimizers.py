"""\
 Minimizers.py: Geometry Minimizers

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Todo list:
# * conjugate gradient geometry minimizer
# * numerical hessian/frequencies

import copy
from math import sqrt
from PyQuante.NumWrap import array

def shift_geo(atomlist,atom,dir,amount):
    atomlist[atom].r[dir] += amount
    return atomlist

def RHFEnergyFunction(atomlist,**kwargs):
    from hartree_fock import rhf
    return rhf(atomlist,**kwargs)[0]

def NumericEnergyForcesHF(atomlist,**kwargs):
    return RHFEnergyFunction(atomlist,**kwargs),NumericForces(atomlist,RHFEnergyFunction)

def NumericForces(atomlist,EnergyFunction):
    "Return the forces on each atom in atomlist via finite differences"
    Forces = []
    dx = 1e-4
    nat = len(atomlist)
    for i in xrange(nat):
        plus_x_geo = shift_geo(copy.deepcopy(atomlist),i,0,dx)
        minus_x_geo = shift_geo(copy.deepcopy(atomlist),i,0,-dx)
        plus_y_geo = shift_geo(copy.deepcopy(atomlist),i,1,dx)
        minus_y_geo = shift_geo(copy.deepcopy(atomlist),i,1,-dx)
        plus_z_geo = shift_geo(copy.deepcopy(atomlist),i,2,dx)
        minus_z_geo = shift_geo(copy.deepcopy(atomlist),i,2,-dx)
        fx = (EnergyFunction(plus_x_geo)-EnergyFunction(minus_x_geo))/(2*dx)
        fy = (EnergyFunction(plus_y_geo)-EnergyFunction(minus_y_geo))/(2*dx)
        fz = (EnergyFunction(plus_z_geo)-EnergyFunction(minus_z_geo))/(2*dx)
        Forces.append((fx,fy,fz))
    return array(Forces)

def Frms(F):
    "Compute the RMS value of a vector of forces"
    sqsum = 0
    for fx,fy,fz in F: sqsum += fx*fx+fy*fy+fz*fz
    return sqrt(sqsum/len(F))

def SteepestDescent(atomlist,EnergyForces=NumericEnergyForcesHF):
    "Called with a pointer to an energy/force evaluator"
    step = 0.1 # this seems awfully small
    Eold = None
    for i in xrange(50):
        E,F = EnergyForces(atomlist)
        print i,E,Frms(F),step
        for j in xrange(len(atomlist)):
            atomlist[j].r -= step*F[j]
        if Eold:
            dE = E-Eold
            if abs(dE) < 0.001: break
            if dE > 0:
                step *= 0.5
            else:
                step *= 1.2
        Eold = E
    return atomlist

def test():
    from PyQuante import Molecule
    h2 = Molecule('H2',
                  [(1,  (0.00000000,     0.00000000,     0.5)),
                   (1,  (0.00000000,     0.00000000,    -0.5))],
                  units='Angstrom')
    h2opt = SteepestDescent(h2)
    print h2opt
    return

if __name__ == '__main__': test()

    
