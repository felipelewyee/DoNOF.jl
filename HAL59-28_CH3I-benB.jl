using DoNOF

mol = """
0 1
I    -0.20532   -3.36852    1.88377
C     0.06920   -5.47193    1.87275
H     0.07360   -5.80564    0.83983
H     1.01805   -5.68986    2.35292
H    -0.75287   -5.92202    2.42041
"""

bset,p = DoNOF.molecule(mol,"def2-qzvp",spherical=true)

p.title = "HAL59-28_CH3I-benB"

p.ipnof = 9

p.RI = true
p.maxit = 40

p.maxloop = 10

DoNOF.set_ncwo(p,1)

C = DoNOF.read_C(title=p.title)
n = DoNOF.read_n(title=p.title)

DoNOF.energy(bset,p,C=C,n=n,do_hfidr=false,do_m_diagnostic=true)
