using DoNOF
using Lints

Z = [8,1,1]
xyz = [[0.0,0.0,0.116],[0.0,0.749,-0.453],[0.0,-0.749,-0.453]]

@lints mol = Lints.Molecule(Z,xyz)

p = DoNOF.Param(mol,"6-31G",Z,0,1)

p.ista=1

p.RI = false
p.gpu = false
DoNOF.set_ncwo(p,1)
#p.HighSpin = true
#p.MSpin = p.nsoc

p.threshl = 10^-5#4   # Convergencia de los multiplicadores de Lagrange
p.threshe = 10^-6#6   # Convergencia de la energía

E,C,gamma,fmiug0 = DoNOF.compute_energy(mol,"6-31G",Z,xyz,p,do_hfidr=true,do_nofmp2=true)
