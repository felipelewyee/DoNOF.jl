using DoNOF

DoNOF.molecule("""
0 1
O  0.0000   0.000   0.116
H  0.0000   0.749  -0.453
H  0.0000  -0.749  -0.453
""")

bas_name = "6-31g"
p = DoNOF.Param(bas_name)

p.RI = true
p.gpu = true
E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,do_hfidr=true)

p.RI = false
E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,C=C,fmiug0=fmiug0,gamma=gamma,do_hfidr=false)
