using DoNOF
using LinearAlgebra

mol = """
0 1
 O  0.0000   0.000   0.121
 H  0.0000   0.751  -0.485
 H  0.0000  -0.751  -0.485
"""

bset, p = DoNOF.molecule(mol, "cc-pvdz", spherical = true)

p.ipnof = 8

p.RI = true
p.gpu = false

p.orb_method = "Rotations"

S, T, V, H, I, b_mnl = DoNOF.compute_integrals(bset, p)

Ei, C = eigen(H, S)

gamma = DoNOF.guess_gamma_trigonometric(p)

n, dn_dgamma = DoNOF.ocupacion_trigonometric(
    gamma,
    p.no1,
    p.ndoc,
    p.nalpha,
    p.nv,
    p.nbf5,
    p.ndns,
    p.ncwo,
    p.HighSpin,
)
cj12, ck12 = DoNOF.PNOFi_selector(n, p)
