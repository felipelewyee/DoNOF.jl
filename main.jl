using DoNOF

mol = """
  O    0.0000000    0.0099701   -0.0000000
  H   -0.0000000   -0.3962400    0.7727381
  H    0.0000000   -0.3962400   -0.7727381
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz")

p.ipnof = 7
p.ista = 1

p.RI = true
p.gpu = true

p.method = "ID"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_m_diagnostic=true)
