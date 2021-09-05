using DoNOF

DoNOF.molecule("""
O  0.0000   0.000   0.116
H  0.0000   0.749  -0.453
H  0.0000  -0.749  -0.453
""")

bas_name = "cc-pvdz"
p = DoNOF.Param(bas_name,0,1)

p.RI = true
p.gpu = true
DoNOF.set_ncwo(p,3)

E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,do_hfidr=true)

p.RI = false
p.gpu = false
E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,C=C,fmiug0=fmiug0,gamma=gamma,do_hfidr=false)
