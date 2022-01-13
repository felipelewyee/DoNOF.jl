using DoNOF

mol = """
O  0.0000   0.000   0.116
H  0.0000   0.749  -0.453
H  0.0000  -0.749  -0.453
"""

bset,p = DoNOF.molecule(mol,"cc-pvtz")

p.gradient = "analytical"

p.ipnof = 7

p.RI = true
p.gpu = true
E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_ekt=true,do_mulliken_pop=true,do_lowdin_pop=true,do_m_diagnostic=true)

#p.RI = false
#E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,C=C,fmiug0=fmiug0,gamma=gamma,do_hfidr=false)
