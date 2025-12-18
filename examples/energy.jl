using DoNOF

mol = """
0 1
 H   0.0 0.0  -0.3804
 H   0.0 0.0   0.3804
"""

bset, p = DoNOF.molecule(mol, "cc-pvdz", spherical = true)

p.RI = true

p.ipnof = 5
p.ista = 0

p.threshgorb = 10^-5

DoNOF.energy(bset, p)
