using Tullio, TensorOperations, LoopVectorization
using BenchmarkTools

N = 200

b = rand(N, N, 7 * N);
C = rand(N, N);

println("RI First MO Transformation")
@benchmark @tullio b_new[i, n, k] := $b[m, n, k] * $C[m, i]
@benchmark @tensor b_new[i, n, k] := $b[m, n, k] * $C[m, i]

println("RI Second MO Transformation")
@benchmark @tullio b_new[m, j, k] := $b[m, n, k] * $C[n, j]
@benchmark @tensor b_new[m, j, k] := $b[m, n, k] * $C[n, j]

println("Density MO Transformation")
@benchmark @tullio D[i, m, n] := C[m, i] * C[n, i]

b_q = rand(N, 7 * N);
println("J for Orbitals")
@benchmark @tullio J[q, m, n] := $b_q[q, l] * $b[m, n, l]
@benchmark @tensor J[q, m, n] := $b_q[q, l] * $b[m, n, l]

println("K for Orbitals")
@benchmark @tullio K[q, m, n] := $b[q, m, l] * $b[q, n, l]
#@benchmark @tensor K[q,m,n] := $b[q,m,l]*$b[q,n,l] This can not be done directly with tensor

println("J for Occupations")
@benchmark @tullio J_MO[p, q] := $b[p, p, l] * $b[q, q, l]
#@benchmark @tensor J_MO[p,q] := $b[p,p,l]*$b[q,q,l] This can not be done directly with tensor

println("K for Occupations")
@benchmark @tullio K_MO[p, q] := $b[p, q, l] * $b[p, q, l]
#@benchmark @tensor K_MO[p,q] := $b[p,q,l]*$b[p,q,l] This can not be done directly with tensor

#QHMATm
D = rand(N, N, N);
H = rand(N, N);
println("H for Occupations")
@benchmark @tullio H_core[i] := $D[i, m, n] * $H[m, n]
@benchmark @tensor H_core[i] := $D[i, m, n] * $H[m, n]
