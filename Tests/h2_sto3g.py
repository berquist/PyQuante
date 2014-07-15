#!/usr/bin/env python
"H2 using Gaussians"

from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule

energy = -1.082098
name = "H2"

def main():
    h2 = Molecule('h2',atomlist=[(1,(0,0,0.6921756113231793)),(1,(0,0,-0.6921756113231793))])
    en,orbe,orbs = rhf(h2,basis="sto3g")
    return en

if __name__ == '__main__':
    print main()

