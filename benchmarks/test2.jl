using DoNOF
using LinearAlgebra

mol = """
0 1
  C    0.0000000    0.6985454   -3.6450775
  C    0.0000000    1.4019537   -2.4373222
  C    0.0000000    0.7082257   -1.2166184
  C   -0.0000000   -0.7082257   -1.2166184
  C   -0.0000000   -1.4019537   -2.4373222
  C    0.0000000   -0.6985454   -3.6450775
  C    0.0000000    1.4063374   -0.0000000
  C    0.0000000    0.7082257    1.2166184
  C   -0.0000000   -0.7082257    1.2166184
  C   -0.0000000   -1.4063374    0.0000000
  H    0.0000000    1.2377776   -4.5835677
  H    0.0000000    2.4853073   -2.4576696
  H   -0.0000000   -2.4853073   -2.4576696
  H    0.0000000   -1.2377776   -4.5835677
  H    0.0000000    2.4911563   -0.0000000
  C    0.0000000    1.4019537    2.4373222
  C   -0.0000000   -1.4019537    2.4373222
  H   -0.0000000   -2.4911563    0.0000000
  C   -0.0000000    0.6985454    3.6450775
  C   -0.0000000   -0.6985454    3.6450775
  H    0.0000000    2.4853073    2.4576696
  H   -0.0000000   -2.4853073    2.4576696
  H   -0.0000000    1.2377776    4.5835677
  H   -0.0000000   -1.2377776    4.5835677
  """

bset, p = DoNOF.molecule(mol, "cc-pvtz", spherical = true)

p.ipnof = 8

p.RI = true

p.occ_method = "Softmax"

S, T, V, H, I, b_mnl = DoNOF.compute_integrals(bset, p)

Ei, C = eigen(H, S)

p.nv = p.nbf5 - p.no1 - p.nsoc
gamma = DoNOF.guess_gamma_softmax(p.ndoc, p.ncwo)

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

@benchmark DoNOF.orbopt_adabelief(gamma, C, H, I, b_mnl, p)
