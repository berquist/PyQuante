from PyQuante.TestMolecules import h2
from PyQuante.PyQuante2 import SCF

hf = SCF(h2,method="HF",basis="3-21g")
print h2
print hf.basis_set.bfs.bfs
print hf.DoAveraging
hf.iterate()
print hf.iterator
print hf.energy, "(benchmark = -1.122956)"
