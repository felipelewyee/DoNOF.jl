function compute_integrals(bset, p)

    # Overlap, Kinetics, Potential
    S = overlap(bset)
    S = (S .+ S') ./ 2
    T = kinetic(bset)
    T = (T .+ T') ./ 2
    V = nuclear(bset)
    V = (V .+ V') ./ 2
    H = T + V
    I = nothing
    println("Basis Set                                            = ", bset.name)
    if (!p.RI)
        # Integrales de Repulsión Electrónica, ERIs (mu nu | sigma lambda)
        I = ERI_2e4c(bset)
    else
        aux = nothing
        if p.spherical
            try
                aux = BasisSet(bset.name * "-jkfit", bset.atoms)
            catch
                aux = BasisSet("def2-universal-jkfit", bset.atoms)
            end
        else
            try
                aux = BasisSet(
                    bset.name * "-jkfit",
                    bset.atoms,
                    spherical = false,
                    lib = :acsint,
                )
            catch
                aux = BasisSet(
                    "def2-universal-jkfit",
                    bset.atoms,
                    spherical = false,
                    lib = :acsint,
                )
            end
        end
        println("Auxiliary Basis Set                                  = ", aux.name)
        println(
            "Gaussian Type                                        = ",
            p.spherical ? "Spherical" : "Cartesian",
        )

        G = ERI_2e2c(aux)
        G = (G .+ G') ./ 2
        I = ERI_2e3c(bset, aux)
        ext = Base.get_extension(@__MODULE__, :DoNOFGPUExt)
        if !isnothing(ext)
            I = ext.Iaux_gpu(G, I)
        else
            I = Iaux(G, I)
        end
        p.nbfaux = size(I)[3]
    end

    ext = Base.get_extension(@__MODULE__, :DoNOFGPUExt)
    if !isnothing(ext)
        I = ext.eris_to_gpu(I)
    end

    return S, T, V, H, I

end

function Iaux(G::Matrix{Float64}, I::Array{Float64,3})

    nbf = size(I)[1]
    nbfaux = size(I)[3]

    G = Symmetric(G)
    evals, evecs = eigen(G)
    sqrtinv = Float64[]
    for i = 1:nbfaux
        if (evals[i] < 0.0)
            append!(sqrtinv, 0.0)
        else
            append!(sqrtinv, 1 / sqrt(evals[i]))
        end
    end
    Gmsqrt = evecs * Diagonal(sqrtinv) * evecs'

    #@tullio I[m, n, l] := I[m, n, k] * Gmsqrt[k, l]
    Threads.@threads for m = 1:nbf
        I[m, :, :] = I[m, :, :] * Gmsqrt
    end

    return I
end

######################################### J_mn^(j) K_mn^(j) #########################################

function computeJKj(C, I::Array{Float64,4}, p)

    nbf5 = p.nbf5

    Cnbf5 = view(C, :, 1:nbf5)

    @tullio D[i, m, n] := Cnbf5[m, i] * Cnbf5[n, i]
    @tullio J[i, m, n] := D[i, s, l] * I[m, n, s, l]
    @tullio K[i, m, s] := D[i, n, l] * I[m, n, s, l]

    return J, K

end

function computeJKj(C, b_mnl::Array{Float64,3}, p)

    nbf5 = p.nbf5

    Cnbf5 = view(C, :, 1:nbf5)

    #b_transform
    @tullio b_qnl[q, n, l] := Cnbf5[m, q] * b_mnl[m, n, l]
    @tullio b_qql[q, l] := Cnbf5[n, q] * b_qnl[q, n, l]

    #hstarj
    @tullio J[q, m, n] := b_qql[q, l] * b_mnl[m, n, l]

    #hstark
    @tullio K[q, m, n] := b_qnl[q, m, l] * b_qnl[q, n, l]

    return J, K

end

######################################### J_pq K_pq #########################################

function computeJKH_MO(C, H, I, p)

    if (p.RI)
        J_MO, K_MO, H_core = JKH_MO(C, H, I, p.nbf, p.nbf5, p.nbfaux)
    else
        J_MO, K_MO, H_core = JKH_MO(C, H, I, p.nbf, p.nbf5)
    end
    return J_MO, K_MO, H_core

end

#########################################

function JKH_MO(
    C::Matrix{Float64},
    H::Matrix{Float64},
    b_mnl::Array{Float64,3},
    nbf::Int64,
    nbf5::Int64,
    nbfaux::Int64,
)

    Cnbf5 = view(C, :, 1:nbf5)

    #denmatj
    @tullio grad=false D[i, m, n] := Cnbf5[m, i] * Cnbf5[n, i]

    #b transform
    @tullio grad=false b_pnl[p, n, l] := Cnbf5[m, p] * b_mnl[m, n, l]
    @tullio grad=false b_pql[p, q, l] := Cnbf5[n, q] * b_pnl[p, n, l]

    #QJMATm
    @tullio grad=false J_MO[p, q] := b_pql[p, p, l] * b_pql[q, q, l]

    #QKMATm
    @tullio grad=false K_MO[p, q] := b_pql[p, q, l] * b_pql[p, q, l]

    #QHMATm
    @tullio grad=false H_core[i] := D[i, m, n] * H[m, n]

    return J_MO, K_MO, H_core
end

function JKH_MO(
    C::Matrix{Float64},
    H::Matrix{Float64},
    I::Array{Float64,4},
    nbf::Int64,
    nbf5::Int64,
)


    Cnbf5 = view(C, :, 1:nbf5)

    #denmatj
    @tullio grad=false D[i, m, n] := Cnbf5[m, i] * Cnbf5[n, i]

    #QJMATm
    @tullio grad=false J[j, m, n] := D[j, s, l] * I[m, n, s, l]
    @tullio grad=false J_MO[i, j] := D[i, m, n] * J[j, m, n]

    #QKMATm
    @tullio grad=false K[j, m, s] := D[j, n, l] * I[m, n, s, l]
    @tullio grad=false K_MO[i, j] := D[i, m, s] * K[j, m, s]

    #QHMATm
    @tullio grad=false H_core[i] := D[i, m, n] * H[m, n]

    return J_MO, K_MO, H_core
end

function computeD_HF(C, I, b_mnl, p)

    Cbeta = view(C, :, 1:p.nbeta)

    @tullio D[m, n] := Cbeta[m, j] * Cbeta[n, j]

    return D
end

function computeDalpha_HF(C, I, b_mnl, p)

    Calpha = view(C, :, (p.nbeta+1):p.nalpha)
    @tullio D[m, n] := Calpha[m, j] * Calpha[n, j]

    return D
end

function computeJK_HF(D, I, p)

    #if(p.gpu)
    #    J,K = JK_HF_Full(CuArray(D),I,p)
    #    return Array(J),Array(K)
    #else
    J, K = JK_HF_Full(D, I, p)
    return J, K
    #end


end

function JK_HF_Full(D, I, p)

    #denmatj
    @tullio J[m, n] := D[l, s] * I[m, n, s, l]
    @tullio K[m, s] := D[n, l] * I[m, n, s, l]

    return J, K

end

function compute_iajb(C, I, p)

    #if(p.gpu)
    #    iajb = iajb_Full(CuArray(C),I,p.no1,p.nalpha,p.nbf,p.nbf5)
    #	return Array(iajb)
    #else
    iajb = iajb_Full(C, I, p.no1, p.nalpha, p.nbf, p.nbf5)
    return iajb
    #end

end

function iajb_Full(C, I, no1, nalpha, nbf, nbf5)

    Cocc = view(C, :, (no1+1):nalpha)
    Ccwo = view(C, :, (nalpha+1):nbf)

    @tullio iajb[i, a, j, b] :=
        ((Ccwo[n, a] * ((Cocc[m, i] * I[m, n, s, l]) * Cocc[s, j])) * Ccwo[l, b])

    return iajb

end

function compute_eris_full(I_AO::Array{Float64,4}, C, nbf5)

    C_nbf5 = @view C[1:end, 1:nbf5]

    @tullio Iinsl[i, n, s, l] := I_AO[m, n, s, l] * C_nbf5[m, i]
    @tullio Iijsl[i, j, s, l] := Iinsl[i, n, s, l] * C_nbf5[n, j]
    Iinsl = nothing
    @tullio Iijkl[i, j, k, l] := Iijsl[i, j, s, l] * C_nbf5[s, k]
    Iijsl = nothing
    @tullio I_MO[i, j, k, r] := Iijkl[i, j, k, l] * C_nbf5[l, r]
    Iijkl = nothing
    #if(pp.gpu):
    #    I = I.get()

    return I_MO
end

function compute_eris_full(b_mnl::Array{Float64,3}, C, nbf5)

    C_nbf5 = @view C[1:end, 1:nbf5]

    @tullio b_pnl[p, n, l] := C_nbf5[m, p] * b_mnl[m, n, l]
    @tullio b_pql[p, q, l] := C_nbf5[n, q] * b_pnl[p, n, l]
    b_pnl = nothing
    @tullio I_MO[p, q, s, l] := b_pql[p, q, R] * b_pql[s, l, R]
    b_pql = nothing

    return I_MO
end
