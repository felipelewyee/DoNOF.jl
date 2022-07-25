using DoNOF

mol = """
  F    0.0000000    0.0000000    0.4927565
  H    0.0000000    0.0000000   -0.4927565
"""

bset,p = DoNOF.molecule(mol,"cc-pvtz")

p.gradient = "analytical"

p.ipnof = 7

p.RI = true
p.gpu = true

p.method = "Rotations"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_ekt=true,do_mulliken_pop=true,do_lowdin_pop=true,do_m_diagnostic=true)
