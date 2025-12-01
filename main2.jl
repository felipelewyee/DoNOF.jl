using DoNOF

mol = """
0 1
 H  0.0000   0.000  -0.3805
 H  0.0000   0.000   0.3805
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz",spherical=true)

p.RI = true

p.ipnof = 8
p.ista = 0

p.threshgorb = 10^-6
#p.threshgocc = 10^-6

DoNOF.energy(bset,p,do_gradients=true)
