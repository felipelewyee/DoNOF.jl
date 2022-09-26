using DoNOF

mol = """
  C   -6.4216976    1.2765240    0.0000017
  C   -5.6078111    2.5592807    0.0000016
  H   -5.9185483    0.5113038   -0.6278778
  H   -7.4350360    1.4734886   -0.4090266
  H   -6.5167561    0.8908563    1.0369098
  H   -5.5127526    2.9449483   -1.0369064
  H   -4.5944727    2.3623161    0.4090299
  H   -6.1109604    3.3245009    0.6278811
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz",spherical=true)

p.ipnof = 8

p.RI = true
p.gpu = true

p.method = "ID"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_m_diagnostic=true)
