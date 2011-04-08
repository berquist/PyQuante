from PyQuante.TestMolecules import h2
from PyQuante.PyQuante2 import SCF

hf = SCF(h2,method="HF",basis="3-21g")
hf.iterate()
print hf.energy, "(benchmark = -1.122956)"
