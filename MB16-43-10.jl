using DoNOF

mol = """
0 2
H     1.88949   -0.89321    1.47608
Si   -1.10595    1.72151   -0.22304
H     1.12997   -0.90678   -1.31889
Cl   -1.30646   -2.33408   -0.31087
C     2.04583   -0.26887    0.58701
H    -2.41985   -0.81980   -1.39550
F    -2.64560    0.05287   -1.69479
H     0.86583    1.41180   -1.67699
C    -1.20888    0.22645    1.08182
N     0.41533    1.04085   -0.84097
B    -0.48686   -0.82592    0.24604
H    -0.57135    0.62992    1.88853
Mg    1.73214    1.83957    0.97924
C     0.87273   -0.35310   -0.40909
H     3.00046   -0.55679    0.13379
H    -2.20684    0.03559    1.47763
"""

bset,p = DoNOF.molecule(mol,"def2-qzvp",spherical=true)

p.title = "MB16-43-10"

p.ipnof = 8

p.RI = true
p.gpu = false

p.maxloop = 10

DoNOF.set_ncwo(p,1)

C = DoNOF.read_C(title=p.title)
n = DoNOF.read_n(title=p.title)

DoNOF.energy(bset,p,C=C,n=n,do_hfidr=false,do_m_diagnostic=true)
