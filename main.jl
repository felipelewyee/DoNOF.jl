using DoNOF

mol = """
0 1
 O  0.0000   0.000   0.121
 H  0.0000   0.751  -0.485
 H  0.0000  -0.751  -0.485
"""

bset,p = DoNOF.molecule(mol,"def2-tzvpd",spherical=true)

p.ipnof = 7

p.RI = true
p.gpu = false

p.orb_method = "ADABelief"
p.maxit = 100

DoNOF.energy(bset,p,do_hfidr=false)
