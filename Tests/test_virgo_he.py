from PyQuante import Molecule, SCF

# Define the atom
he = Molecule('he',atomlist=[(2,(0,0,0))],charge=0,multiplicity=1)

# Run HF
solver = SCF(he,method="HF",basis="dzvp")
solver.iterate()

# Show result
print solver
print "HF Energy = ",solver.energy

print "---With fractional charge---"
he = Molecule('he',atomlist=[(2.04,(0,0,0))],charge=0,multiplicity=1)

# Run HF
solver = SCF(he,method="HF",basis="dzvp")
solver.iterate()

# Show result
print solver
print "HF Energy = ",solver.energy
