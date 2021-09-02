using DoNOF

DoNOF.molecule("""
O 0.0 0.0 0.116
H 0.0 0.749 -0.453
H 0.0 -0.749 -0.453
""")

DoNOF.molecule_from_file("h2o.xyz")

bas_name = "cc-pvtz"
p = DoNOF.Param(bas_name,0,1)

#p.ista=1

p.RI = true
p.gpu = false#true
DoNOF.set_ncwo(p,1)
#p.HighSpin = true
#p.MSpin = p.nsoc

p.threshl = 10^-5#4   # Convergencia de los multiplicadores de Lagrange
p.threshe = 10^-6#6   # Convergencia de la energía

E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,do_hfidr=true)

p.RI = false
p.gpu = true
E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,C=C,fmiug0=fmiug0,gamma=gamma,do_hfidr=false,do_nofmp2=true,nofmp2strategy="numerical")
