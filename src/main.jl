using DoNOF

mol = """
O  0.0000   0.000   0.116
H  0.0000   0.749  -0.453
H  0.0000  -0.749  -0.453
"""

bset,p = DoNOF.molecule(mol,"cc-pvtz")

p.gradient = "analytical"

#p.ipnof = 7

p.RI = false
p.gpu = true
p.EBI = true
DoNOF.set_muller(p,true)

p.method = "Rotations"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=false)
#,do_ekt=true,do_mulliken_pop=true,do_lowdin_pop=true,do_m_diagnostic=true)
