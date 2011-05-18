from PyQuante.TestMolecules import h2
from PyQuante.LA2 import geigh,mkdens
from PyQuante.Ints import get2JmK,getbasis,getints
from PyQuante.Convergence import DIIS
from PyQuante.hartree_fock import get_energy

bfs = getbasis(h2,basis="3-21g")
nclosed,nopen = h2.get_closedopen()
nocc = nclosed
nel = h2.get_nel()
S,h,Ints = getints(bfs,h2)

orbe,orbs = geigh(h,S)

enuke = h2.get_enuke()
eold = 0

avg = DIIS(S)

for i in range(20):
    D = mkdens(orbs,0,nocc)
    G = get2JmK(Ints,D)
    F = h+G
    F = avg.getF(F,D)
    orbe,orbs = geigh(F,S)
    energy = get_energy(h,F,D,enuke)
    print i,energy,energy-eold
    if abs(energy-eold)<1e-5:
        break
    eold = energy
print "Converged"
print energy, "(benchmark = -1.122956)"
