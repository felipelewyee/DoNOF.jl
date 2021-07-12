using PyCall
using Tullio

include("parameters.jl")
include("energy.jl")

psi4 = pyimport_conda("psi4", "psi4")
np = pyimport_conda("numpy","psi4")

mol = psi4.geometry("""
0 1
O  0.0000   0.000   0.116
H  0.0000   0.749  -0.453
H  0.0000  -0.749  -0.453
  symmetry c1
""")
mol = psi4.geometry("""
0 1
Mn  0.0000   0.000   0.000
Mn  0.0000   0.000   3.140
  symmetry c1
""")

psi4.set_options(Dict("basis"=>"aug-cc-pVDZ"))

wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("basis"))

#p = parameters.Param(mol.natom(),wfn.basisset().nbf(),wfn.nalpha(),wfn.nbeta(),mol.multiplicity(),Z)
p = parameters.Param(mol,wfn)
p.RI = true
p.gpu = true
parameters.set_ncwo(p,1)
#p.HighSpin = true
#p.MSpin = p.nsoc

energy.compute_energy(wfn,mol,p,hfidr=false)
