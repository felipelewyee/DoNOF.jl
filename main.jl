using DoNOF

mol = """
0 1
  O   -0.0000000   -0.0072090    0.0000000
  H   -0.0000000    0.3369355    0.7705103
  H   -0.0000000    0.3369355   -0.7705103
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz",spherical=false)

p.ipnof = 8

p.RI = false
p.gpu = false

p.method = "ID"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_m_diagnostic=true,do_mbpt=false,do_translate_to_donofsw=true)
