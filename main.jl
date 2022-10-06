using DoNOF

mol = """
  O   -0.0000000   -0.0155252    0.0000000
  H    0.0000000    0.4944556   -0.7830366
  H   -0.0000000    0.4944556    0.7830366
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz",spherical=true)

p.ipnof = 8

p.RI = true
p.gpu = true

p.method = "ID"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_m_diagnostic=true)
