"""\
Defaults.py - PyQuante-wide program defaults.
"""

# SCF flags
MaxIter = 30
DerivativeType = 'analyt'
MolecularCharge = 0
SpinMultiplicity = 1
LengthUnits = 'bohr'
ConvergenceCriteria = 1e-5
HamiltonianMethod = 'HF'
Averaging = True
MixingFraction = 0.5
ElectronTemperature = False
FDTolerance = 1e-9
FDOccTolerance = 1e-5

# MINDO flags
MINDOAveraging = False
MINDOReturnEtot = False
MINDONumericForces = False

# DFT-related flags
DFTFunctional = "SVWN"
DFTSpinPolarized = False
DFTConvergenceCriteria = 1e-4
DFTAveraging = Averaging # Consider removing?
DFTElectronTemperature = ElectronTemperature # consider removing?
DFTGridRadii = 32
DFTGridFineness = 1
DFTDensityGradient = False
DFTRadialGridType = 'EulerMaclaurin'
DFTGridSG1 = True
DFTGridAngularPoints = 194
DFTDensityCutoff = 1e-10
AM05DensityCutoff = 1e-16
DFTXalphaFactor = 2./3.
DFTBeckeHetero = True
OEPOptMethod = 'BFGS'
OEPIters = 100
OEPTolerance = 1e-5

# Density matrix purification stuff
DMPTolerance = 1e-7
DMPIterations = 100
DMPLanczosMinmaxIters = 8

# Integral stuff
IntsOmitF = False

# Math stuff
OrthogMethod = 'Chol'
NumericForceDx = 1e-6
NumericForceSym = True
SubspaceVirtualOrbs = 1
DavidsonEvecTolerance = 1e-6
DavidsonNormTolerance = 1e-10
JacobiSweeps = 100
JacobiTolerance = 1e-10
ExpSteps = 12
ExpCutoff = 1e-8

# Dynamics options
DynSteps = 100
DynJob = 'pydyn'
DynTStep = 0.1
