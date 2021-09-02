using DoNOF

DoNOF.molecule("""
  O    0.0000000    0.0184016   -0.0000000
  H    0.0000000   -0.5383158   -0.7588684
  H   -0.0000000   -0.5383158    0.7588684
""")

bas_name = "cc-pvdz"
p = DoNOF.Param(bas_name,0,1)

#p.ista=1

p.RI = true
p.gpu = true
DoNOF.set_ncwo(p,1)
#p.HighSpin = true
#p.MSpin = p.nsoc

p.threshl = 10^-5#4   # Convergencia de los multiplicadores de Lagrange
p.threshe = 10^-6#6   # Convergencia de la energía

E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,do_hfidr=true)

p.RI = false
p.gpu = false
E,C,gamma,fmiug0 = DoNOF.compute_energy(bas_name,p,C=C,fmiug0=fmiug0,gamma=gamma,do_hfidr=false,do_nofmp2=false,nofmp2strategy="numerical")
