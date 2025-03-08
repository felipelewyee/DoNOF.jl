using CUDA, KernelAbstractions, cuTENSOR
using Tullio, TensorOperations, LoopVectorization
using BenchmarkTools

N = 200

b = rand(N, N, 7 * N);
C = rand(N, N);
b = CuArray{Float64}(b);
C = CuArray{Float64}(C);

println("RI First MO Transformation")
CUDA.@time @tullio b_new[i, n, k] := $b[m, n, k] * $C[m, i];
CUDA.@time @tensor b_new[i, n, k] := b[m, n, k] * C[m, i];

println("RI Second MO Transformation")
CUDA.@time @tullio b_new[m, j, k] := $b[m, n, k] * $C[n, j];
CUDA.@time @tensor b_new[m, j, k] := b[m, n, k] * C[n, j];

println("Density MO Transformation")
CUDA.@time @tullio D[i, m, n] := C[m, i] * C[n, i];

b_q = rand(N, 7 * N);
b_q = CuArray{Float64}(b_q);
println("J for Orbitals")
CUDA.@time @tullio J[q, m, n] := $b_q[q, l] * $b[m, n, l];
CUDA.@time @tensor J[q, m, n] := b_q[q, l] * b[m, n, l];

println("K for Orbitals")
CUDA.@time @tullio K[q, m, n] := $b[q, m, l] * $b[q, n, l];
#@benchmark @tensor K[q,m,n] := $b[q,m,l]*$b[q,n,l] This can not be done directly with tensor

println("J for Occupations")
CUDA.@time @tullio J_MO[p, q] := $b[p, p, l] * $b[q, q, l];
#@benchmark @tensor J_MO[p,q] := $b[p,p,l]*$b[q,q,l] This can not be done directly with tensor

println("K for Occupations")
CUDA.@time @tullio K_MO[p, q] := $b[p, q, l] * $b[p, q, l];
#@benchmark @tensor K_MO[p,q] := $b[p,q,l]*$b[p,q,l] This can not be done directly with tensor

#QHMATm
D = rand(N, N, N);
D = CuArray{Float64}(D);
H = rand(N, N);
H = CuArray{Float64}(H);
println("H for Occupations")
CUDA.@time @tullio H_core[i] := $D[i, m, n] * $H[m, n];
CUDA.@time @tensor H_core[i] := D[i, m, n] * H[m, n];
