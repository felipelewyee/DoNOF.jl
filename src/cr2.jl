using DoNOF

mol = """
Cr  0.0000  0.0000  0.0000
Cr  1.6600  0.0000  0.0000
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz")

p.gradient = "analytical"


#p.RI = true
p.gpu = true
p.EBI = true
DoNOF.set_muller(p,true)

p.method = "Rotations"

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=false)
#,do_ekt=true,do_mulliken_pop=true,do_lowdin_pop=true,do_m_diagnostic=true)
