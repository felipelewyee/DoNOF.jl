using BenchmarkTools

nbf = 100
nbfaux = 300
nbf5 = 40
nalpha = 20
nbeta = 20

C = rand(nbf,nbf)
n = rand(nbf5)
H = rand(nbf,nbf)
b_mnl = rand(nbf,nbf,nbfaux)
cj12 = rand(nbf5,nbf5)
ck12 = rand(nbf5,nbf5)

@benchmark DoNOF.compute_Lagrange2($C, $n, $H, $b_mnl, $cj12, $ck12, $nbf5, $nalpha, $nbeta)
