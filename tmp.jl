using BenchmarkTools

nbf = 200
nbfaux = 600
nbf5 = 90
nalpha = 30
nbeta = 30

C = rand(nbf,nbf)
n = rand(nbf5)
H = rand(nbf,nbf)
b_mnl = rand(nbf,nbf,nbfaux)
cj12 = rand(nbf5,nbf5)
ck12 = rand(nbf5,nbf5)

@benchmark DoNOF.compute_Lagrange2($C, $n, $H, $b_mnl, $cj12, $ck12, $nbf5, $nalpha, $nbeta)
