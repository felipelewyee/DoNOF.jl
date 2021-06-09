using PyCall, Conda
using Tullio

include("parameters.jl")
include("energy.jl")

psi4 = pyimport_conda("psi4", "psi4")
np = pyimport_conda("numpy","psi4")

mol = psi4.geometry("""
0 3
O  0.0000   0.000   0.116
H  0.0000   0.749  -0.453
H  0.0000  -0.749  -0.453
  symmetry c1
""")

psi4.set_options(Dict("basis"=>"aug-cc-pVDZ"))

wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("basis"))

Z = [mol.Z(i-1) for i=1:mol.natom()]

p = parameters.Param(mol.natom(),wfn.basisset().nbf(),wfn.nalpha(),wfn.nbeta(),mol.multiplicity(),Z)
parameters.set_ncwo(p,1)

# Integrador
mints = psi4.core.MintsHelper(wfn.basisset())

# Overlap, Kinetics, Potential
S = copy(np.asarray(mints.ao_overlap()))
T = copy(np.asarray(mints.ao_kinetic()))
V = copy(np.asarray(mints.ao_potential()))
H = T + V
I = []
b_mnl = []
if (!p.RI)
    # Integrales de Repulsión Electrónica, ERIs (mu nu | sigma lambda)
    I = copy(np.asarray(mints.ao_eri()))
else

    orb = wfn.basisset()
    aux = psi4.core.BasisSet.build(mol, "DF_BASIS_SCF", "", "JKFIT", orb.blend())
    zero_bas = psi4.core.BasisSet.zero_ao_basis_set()

    Ppq = mints.ao_eri(orb, orb, aux, zero_bas)

    metric = mints.ao_eri(aux, zero_bas, aux, zero_bas)
    metric.power(-0.5, 1.e-14)
    p.nbfaux = metric.shape[0]

    Ppq = np.squeeze(Ppq)
    metric = copy(np.squeeze(metric))

    @tullio b_mnl := Ppq[p,q,P]*metric[P,Q]
end

E_nuc = mol.nuclear_repulsion_energy()

energy.compute_energy(S,T,V,H,I,b_mnl,E_nuc,p,hfidr=true)
