using DoNOF

mol = """
0 1
  O  0.0000   0.000   0.121
  H  0.0000   0.751  -0.485
  H  0.0000  -0.751  -0.485
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz",spherical=true)

p.ipnof = 8

p.RI = true
p.gpu = true

p.method = "ID"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_m_diagnostic=true,do_mbpt=false,do_translate_to_donofsw=true)
