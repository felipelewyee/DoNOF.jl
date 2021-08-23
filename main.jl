using PyCall
using Tullio

using DoNOF

psi4 = pyimport_conda("psi4", "psi4")
np = pyimport_conda("numpy","psi4")

mol = psi4.geometry("""
0 1
O  0.0000   0.000   0.116
H  0.0000   0.749  -0.453
H  0.0000  -0.749  -0.453
  symmetry c1
""")
#mol = psi4.geometry("""
#0 1
#Mn  0.0000   0.000   0.000
#Mn  0.0000   0.000   3.150
#  symmetry c1
#""")

psi4.set_options(Dict("basis"=>"6-31G"))

wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("basis"))

#p = parameters.Param(mol.natom(),wfn.basisset().nbf(),wfn.nalpha(),wfn.nbeta(),mol.multiplicity(),Z)
p = DoNOF.Param(mol,wfn)

p.ista=1

p.RI = false
p.gpu = false
DoNOF.set_ncwo(p,1)
#p.HighSpin = true
#p.MSpin = p.nsoc

p.threshl = 10^-5#4   # Convergencia de los multiplicadores de Lagrange
p.threshe = 10^-6#6   # Convergencia de la energía


E,C,gamma,fmiug0 = DoNOF.compute_energy(wfn,mol,p,do_hfidr=true,do_nofmp2=true)
