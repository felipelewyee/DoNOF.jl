module DoNOFGPUExt

using DoNOF
using TensorOperations
using Tullio
using LinearAlgebra
using CUDA, KernelAbstractions, cuTENSOR

function eris_to_gpu(I)

    device = CUDA.device()
    gpu_name = CUDA.name(device)
    println(gpu_name)

    I_gpu = nothing
    if occursin("TX", gpu_name)
        I_gpu = CuArray{Float32}(I)
    else
        I_gpu = CuArray{Float64}(I)
    end

    return I_gpu
end

function DoNOF.computeJKH_MO(C, H, I::CuArray, p)

    if (p.RI)
        J_MO, K_MO, H_core = DoNOF.JKH_MO(C, H, I, p.nbf, p.nbf5, p.nbfaux)
    else
        J_MO, K_MO, H_core = DoNOF.JKH_MO(C, H, I, p.nbf, p.nbf5)
    end
    return Array(J_MO), Array(K_MO), Array(H_core)

end

function DoNOF.JKH_MO(C, H, b_mnl::CuArray{T, 3}, nbf, nbf5, nbfaux) where {T <: Union{Float32, Float64}}

    C = CuArray{typeof(b_mnl).parameters[1]}(C)
    H = CuArray{typeof(b_mnl).parameters[1]}(H)

    Cnbf5 = C[1:nbf, 1:nbf5]

    #b transform
    @tensor b_pnl[p, n, l] := Cnbf5[m, p] * b_mnl[m, n, l]
    @tensor b_pql[p, q, l] := Cnbf5[n, q] * b_pnl[p, n, l]
    CUDA.unsafe_free!(b_pnl)

    #QJMATm
    @tullio tmp[p, l] := b_pql[p, p, l]
    J_MO = tmp * tmp'
    CUDA.unsafe_free!(tmp)
    #@tullio J_MO[p,q] := b_pql[p,p,l]*b_pql[q,q,l]

    #QKMATm
    @tullio K_MO[p, q] := b_pql[p, q, l] * b_pql[p, q, l]
    #K_MO = dropdims( sum(b_pql .* b_pql, dims=3), dims=3)
    CUDA.unsafe_free!(b_pql)

    #QHMATm
    @tensor tmp[m, i] := H[m, n] * Cnbf5[n, i]
    @tullio H_core[i] := tmp[m, i] * Cnbf5[m, i]
    CUDA.unsafe_free!(tmp)
    #@tullio D[i,m,n] := Cnbf5[m,i]*Cnbf5[n,i]
    #@tensor H_core[i] := D[i,m,n]*H[m,n]

    return J_MO, K_MO, H_core
end

function DoNOF.JKH_MO(C, H, I::CuArray{T, 4}, nbf::Int64, nbf5::Int64)  where {T <: Union{Float32, Float64}}

    C = CuArray{typeof(I).parameters[1]}(C)
    H = CuArray{typeof(I).parameters[1]}(H)

    Cnbf5 = view(C, :, 1:nbf5)

    #denmatj
    @tullio grad=false D[i, m, n] := Cnbf5[m, i] * Cnbf5[n, i]

    #QJMATm
    @tensor J[j, m, n] := D[j, s, l] * I[m, n, s, l]
    @tensor J_MO[i, j] := D[i, m, n] * J[j, m, n]

    #QKMATm
    @tensor K[j, m, s] := D[j, n, l] * I[m, n, s, l]
    @tensor K_MO[i, j] := D[i, m, s] * K[j, m, s]

    #QHMATm
    @tensor H_core[i] := D[i, m, n] * H[m, n]

    return J_MO, K_MO, H_core
end

function DoNOF.JKj_RI(C, b_mnl::CuArray, nbf, nbf5, nbfaux)

    Cnbf5 = C[1:nbf, 1:nbf5]
    Cnbf5 = CuArray{typeof(b_mnl).parameters[1]}(Cnbf5)

    #b_transform
    @tensor b_qnl[q, n, l] := Cnbf5[m, q] * b_mnl[m, n, l]

    #hstark
    #@tullio K[q,m,n] := b_qnl[q,m,l]*b_qnl[q,n,l]
    K = CUDA.zeros(typeof(b_mnl).parameters[1], nbf5, nbf, nbf)
    for q = 1:nbf5
        b_q = b_qnl[q, 1:nbf, 1:nbfaux]
        K[q, 1:nbf, 1:nbf] = b_q * b_q'
        CUDA.unsafe_free!(b_q)
    end

    #b_qql = dropdims( sum(Cnbf5' .* b_qnl, dims=2), dims=2)
    @tullio b_qql[q, l] := Cnbf5[n, q] * b_qnl[q, n, l]
    CUDA.unsafe_free!(b_qnl)
    CUDA.unsafe_free!(Cnbf5)

    #hstarj
    @tensor J[q, m, n] := b_qql[q, l] * b_mnl[m, n, l]
    CUDA.unsafe_free!(b_qql)

    return J, K

end

function DoNOF.compute_Lagrange2(C, n, H, b_mnl::CuArray{T, 3}, cj12, ck12, nalpha, nbeta)  where {T <: Union{Float32, Float64}}

    nbf = size(C)[1]
    nbfaux = size(b_mnl)[3]
    nbf5 = size(n)[1]

    C = CuArray{typeof(b_mnl).parameters[1]}(C)
    H = CuArray{typeof(b_mnl).parameters[1]}(H)

    Cnbf5 = C[1:nbf, 1:nbf5]
    H_mat = C' * (H * Cnbf5)
    #@tullio H_mat[i,j] := Cnew[m,i] * H[m,nn] * Cnbf5[nn,j]
    @tensor tmp[m, q, l] := Cnbf5[nn, q] * b_mnl[m, nn, l]
    @tensor b_MO[p, q, l] := C[m, p] * tmp[m, q, l]
    CUDA.unsafe_free!(tmp)

    elag = CUDA.zeros(typeof(b_mnl).parameters[1], nbf, nbf)
    cj12 = CuArray{typeof(b_mnl).parameters[1]}(cj12)
    ck12 = CuArray{typeof(b_mnl).parameters[1]}(ck12)
    n = CuArray{typeof(b_mnl).parameters[1]}(n)

    n_beta = view(n, 1:nbeta)
    n_alpha = view(n, nalpha+1:nbf5)
    Hmat_nbf5 = view(H_mat, 1:nbf, 1:nbf5)
    grad_nbf5 = view(elag, 1:nbf, 1:nbf5)
    grad_nbeta = view(elag, 1:nbf, 1:nbeta)
    grad_nalpha = view(elag, 1:nbf, nalpha+1:nbf5)

    # 2ndH/dy_ab
    elag[1:nbf, 1:nbf5] += Hmat_nbf5 .* n'
    #@tullio grad_nbf5[a,b]  += n[b]*Hmat_nbf5[a,b]

    # dJ_pp/dy_ab
    if nbeta > 0
        b_nbeta_beta = view(b_MO, 1:nbeta, 1:nbeta, 1:nbfaux)
        b_nbf_beta = view(b_MO, 1:nbf, 1:nbeta, 1:nbfaux)

        @tullio tmp[b, k] := b_nbeta_beta[b, b, k]
        tmp2 = n_beta .* tmp
        CUDA.unsafe_free!(tmp)
        tmp3 = permutedims(b_nbf_beta, (2, 3, 1))
        tmp4 = tmp2 .* tmp3
        CUDA.unsafe_free!(tmp2)
        CUDA.unsafe_free!(tmp3)
        tmp5 = sum(tmp4, dims = 2)
        CUDA.unsafe_free!(tmp4)
        tmp6 = dropdims(tmp5, dims = 2)
        CUDA.unsafe_free!(tmp5)
        elag[1:nbf, 1:nbeta] += tmp6'
        CUDA.unsafe_free!(tmp6)
        #@tullio grad_nbeta[a,b] += n_beta[b]*b_nbf_beta[a,b,k]*b_nbeta_beta[b,b,k]

        b_nalpha_alpha =
            view(b_MO, nalpha+1:nbf5, nalpha+1:nbf5, 1:nbfaux)
        b_nbf_alpha = view(b_MO, 1:nbf, nalpha+1:nbf5, 1:nbfaux)
        @tullio tmp[b, k] := b_nalpha_alpha[b, b, k]
        tmp2 = n_alpha .* tmp
        CUDA.unsafe_free!(tmp)
        tmp3 = permutedims(b_nbf_alpha, (2, 3, 1))
        tmp4 = tmp2 .* tmp3
        CUDA.unsafe_free!(tmp2)
        CUDA.unsafe_free!(tmp3)
        tmp5 = sum(tmp4, dims = 2)
        CUDA.unsafe_free!(tmp4)
        tmp6 = dropdims(tmp5, dims = 2)
        CUDA.unsafe_free!(tmp5)
        elag[1:nbf, nalpha+1:nbf5] += tmp6'
        CUDA.unsafe_free!(tmp6)
        #@tullio grad_nalpha[a,b] += n_alpha[b]*b_nbf_alpha[a,b,k]*b_nalpha_alpha[b,b,k]
    end

    b_nbf_nbf5 = view(b_MO, 1:nbf, 1:nbf5, 1:nbfaux)
    b_nbf5_nbf5 = view(b_MO, 1:nbf5, 1:nbf5, 1:nbfaux)

    # C^J_pq dJ_pq/dy_ab
    #@tullio grad_nbf5[a,b] += b_nbf_nbf5[a,b,k]*cj12[b,q]*b_nbf5_nbf5[q,q,k]
    @tullio tmp[q, k] := b_nbf5_nbf5[q, q, k]
    @tensor tmp2[b, k] := cj12[b, q] * tmp[q, k]
    CUDA.unsafe_free!(tmp)
    @tullio grad_nbf5[a, b] += b_nbf_nbf5[a, b, k] * tmp2[b, k]
    CUDA.unsafe_free!(tmp2)

    # -C^K_pq dK_pq/dy_ab
    #@tullio grad_nbf5[a,b] += -ck12[b,q]*b_nbf_nbf5[a,q,k]*b_nbf5_nbf5[b,q,k]
    tmp = ck12 .* b_nbf5_nbf5
    @tensor grad_nbf5[a, b] += -b_nbf_nbf5[a, q, k] * tmp[b, q, k]
    CUDA.unsafe_free!(tmp)
    
    CUDA.unsafe_free!(b_MO)

    elag = Array(elag)
    H_mat = Array(H_mat)

    return elag, H_mat

end


function DoNOF.compute_Lagrange2(C, n, H, I::CuArray{T, 4}, cj12, ck12, nalpha, nbeta)  where {T <: Union{Float32, Float64}}

    nbf = size(C)[1]
    nbf5 = size(n)[1]
    C = CuArray{typeof(I).parameters[1]}(C)
    H = CuArray{typeof(I).parameters[1]}(H)

    Cnbf5 = C[1:nbf, 1:nbf5]
    H_mat = C' * (H * Cnbf5)
    #@tullio H_mat[i,j] := Cnew[m,i] * H[m,nn] * Cnbf5[nn,j]
    @tullio tmp[m, q, s, l] := Cnbf5[nn, q] * I[m, nn, s, l]
    @tullio tmp2[m, q, r, l] := Cnbf5[s, r] * tmp[m, q, s, l]
    @tullio tmp[m, q, r, t] := Cnbf5[l, t] * tmp2[m, q, r, l]
    @tullio I_MO[p, q, r, t] := C[m, p] * tmp[m, q, r, t]

    elag = CUDA.zeros(typeof(I).parameters[1], nbf, nbf)
    cj12 = CuArray{typeof(I).parameters[1]}(cj12)
    ck12 = CuArray{typeof(I).parameters[1]}(ck12)
    n = CuArray{typeof(I).parameters[1]}(n)

    n_beta = view(n, 1:nbeta)
    n_alpha = view(n, nalpha+1:nbf5)
    Hmat_nbf5 = view(H_mat, 1:nbf, 1:nbf5)
    grad_nbf5 = view(elag, 1:nbf, 1:nbf5)
    grad_nbeta = view(elag, 1:nbf, 1:nbeta)
    grad_nalpha = view(elag, 1:nbf, nalpha+1:nbf5)
    I_nb_nb_nb = view(I_MO, 1:nbf, 1:nbeta, 1:nbeta, 1:nbeta)
    I_na_na_na = view(
        I_MO,
        1:nbf,
        nalpha+1:nbf5,
        nalpha+1:nbf5,
        nalpha+1:nbf5,
    )
    I_nbf5_nbf5_nbf5 = view(I_MO, 1:nbf, 1:nbf5, 1:nbf5, 1:nbf5)

    # 2ndH/dy_ab
    @tullio grad_nbf5[a, b] += n[b] * Hmat_nbf5[a, b]

    # dJ_pp/dy_ab
    if nbeta > 0
        @tullio grad_nbeta[a, b] += n_beta[b] * I_nb_nb_nb[a, b, b, b]
        @tullio grad_nalpha[a, b] += n_alpha[b] * I_na_na_na[a, b, b, b]
    end

    # C^J_pq dJ_pq/dy_ab
    @tullio grad_nbf5[a, b] += cj12[b, q] * I_nbf5_nbf5_nbf5[a, b, q, q]

    # -C^K_pq dK_pq/dy_ab
    @tullio grad_nbf5[a, b] += -ck12[b, q] * I_nbf5_nbf5_nbf5[a, q, b, q]
    CUDA.unsafe_free!(I_MO)

    elag = Array(elag)
    H_mat = Array(H_mat)

    return elag, H_mat

end


function DoNOF.computeF_RC(J::CuArray, K, n, H, cj12, ck12, p)

    ini = 0
    if (p.no1 > 1)
        ini = p.no1
    end

    # Matriz de Fock Generalizada
    n = CuArray{typeof(J).parameters[1]}(n)
    H = CuArray{typeof(J).parameters[1]}(H)

    # nH
    @tullio F[i, m, s] := n[i] * H[m, s]

    # -C^K K
    ck12_ini_nbf5 = view(ck12, ini+1:p.nbf5, ini+1:p.nbf5)
    ck12_ini_nbf5[diagind(ck12_ini_nbf5)] .= 0.0
    ck12_gpu = CuArray{typeof(J).parameters[1]}(ck12)
    @tensor F[i, m, n] += -ck12_gpu[i, j] * K[j, m, n]
    CUDA.unsafe_free!(ck12_gpu)
    CUDA.unsafe_free!(K)

    # C^J J
    cj12_ini_nbf5 = view(cj12, ini+1:p.nbf5, ini+1:p.nbf5)
    cj12_ini_nbf5[diagind(cj12_ini_nbf5)] .= 0.0
    cj12_gpu = CuArray{typeof(J).parameters[1]}(cj12)
    @tensor F[i, m, n] += cj12_gpu[i, j] * J[j, m, n]
    CUDA.unsafe_free!(cj12_gpu)

    # nJ
    n_ini_beta = n[ini+1:p.nbeta]
    J_ini_beta = J[ini+1:p.nbeta, 1:p.nbf, 1:p.nbf]
    F[ini+1:p.nbeta, 1:p.nbf, 1:p.nbf] += n_ini_beta .* J_ini_beta
    CUDA.unsafe_free!(n_ini_beta)
    CUDA.unsafe_free!(J_ini_beta)

    n_alpha_nbf5 = n[p.nalpha+1:p.nbf5]
    J_alpha_nbf5 = J[p.nalpha+1:p.nbf5, 1:p.nbf, 1:p.nbf]
    F[p.nalpha+1:p.nbf5, 1:p.nbf, 1:p.nbf] += n_alpha_nbf5 .* J_alpha_nbf5
    CUDA.unsafe_free!(n_alpha_nbf5)
    CUDA.unsafe_free!(J_alpha_nbf5)
    CUDA.unsafe_free!(J)

    return Array(F)

end

function DoNOF.rotate_orbital(y::Array, C, p)

    ynew = zeros(p.nbf, p.nbf)

    n = 1
    for i = 1:p.nbf5
        for j = i+1:p.nbf
            ynew[i, j] = y[n]
            ynew[j, i] = -y[n]
            n += 1
        end
    end

    #tmp = [i<j&&i<=p.nbf5 ? y[Int64((2*p.nbf*i - i^2 - i)/2 + j - p.nbf)] : 0 for i in 1:p.nbf, j in 1:p.nbf]
    #ynew = tmp .- tmp'

    ynew = CuArray{Float32}(ynew)

    vals, vecs = eigen(1im .* ynew)
    evals = exp.(vals ./ 1im)
    U = real.(vecs * Diagonal(evals) * vecs')

    cC = CuArray{Float32}(C)
    @tensor Cnew[m, p] := cC[m, r] * U[r, p]

    return Array(Cnew)

end

end
