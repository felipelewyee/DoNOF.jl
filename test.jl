using DoNOF

DoNOF.molecule("""
0 1
O  0.0000   0.000   0.116
H  0.0000   0.749  -0.453
H  0.0000  -0.749  -0.453
""")

bas_name = "cc-pvdz"
p = DoNOF.Param(bas_name)

p.RI = true
p.gpu = true
E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,do_hfidr=true)

p.RI = false
p.gpu = true
E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,C=C,fmiug0=fmiug0,gamma=gamma,do_hfidr=false,do_nofmp2=true,nofmp2strategy="analytical",tolnofmp2=1e-4)

