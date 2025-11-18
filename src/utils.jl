function compute_Lagrange2(
    C::Matrix{Float64},
    n::Vector{Float64},
    H::Matrix{Float64},
    b_mnl::Array{Float64,3},
    cj12::Matrix{Float64},
    ck12::Matrix{Float64},
    nalpha::Int64,
    nbeta::Int64,
)

    nbf = size(C)[1]
    nbfaux = size(b_mnl)[3]
    nbf5 = size(n)[1]

    Cnbf5 = view(C, 1:nbf, 1:nbf5)
    @tullio grad=false tmp[m, j] := H[m, nn] * Cnbf5[nn, j]
    @tullio grad=false H_mat[i, j] := C[m, i] * tmp[m, j]
    @tullio grad=false tmp[m, q, l] := Cnbf5[nn, q] * b_mnl[m, nn, l]
    @tullio grad=false b_MO[p, q, l] := C[m, p] * tmp[m, q, l]

    elag = zeros(nbf, nbf)

    n_beta = view(n, 1:nbeta)
    n_alpha = view(n, (nalpha+1):nbf5)
    Hmat_nbf5 = view(H_mat, 1:nbf, 1:nbf5)
    grad_nbf5 = view(elag, 1:nbf, 1:nbf5)
    grad_nbeta = view(elag, 1:nbf, 1:nbeta)
    grad_nalpha = view(elag, 1:nbf, (nalpha+1):nbf5)

    b_nbf_beta = view(b_MO, 1:nbf, 1:nbeta, 1:nbfaux)
    b_nbf_alpha = view(b_MO, 1:nbf, (nalpha+1):nbf5, 1:nbfaux)
    b_nbeta_beta = view(b_MO, 1:nbeta, 1:nbeta, 1:nbfaux)
    b_nalpha_alpha = view(b_MO, (nalpha+1):nbf5, (nalpha+1):nbf5, 1:nbfaux)
    b_nbf_nbf5 = view(b_MO, 1:nbf, 1:nbf5, 1:nbfaux)
    b_nbf5_nbf5 = view(b_MO, 1:nbf5, 1:nbf5, 1:nbfaux)

    # 2ndH/dy_ab
    @tullio grad=false grad_nbf5[a, b] += n[b] * Hmat_nbf5[a, b]

    # dJ_pp/dy_ab
    if nbeta > 0
        @tullio grad=false grad_nbeta[a, b] +=
            n_beta[b] * b_nbf_beta[a, b, k] * b_nbeta_beta[b, b, k]
        @tullio grad=false grad_nalpha[a, b] +=
            n_alpha[b] * b_nbf_alpha[a, b, k] * b_nalpha_alpha[b, b, k]
    end

    # C^J_pq dJ_pq/dy_ab 
    @tullio grad=false tmp[b, k] := cj12[b, q] * b_nbf5_nbf5[q, q, k]
    @tullio grad=false grad_nbf5[a, b] += b_nbf_nbf5[a, b, k] * tmp[b, k]

    # -C^K_pq dK_pq/dy_ab 
    @tullio grad=false grad_nbf5[a, b] +=
        -ck12[b, q] * b_nbf_nbf5[a, q, k] * b_nbf5_nbf5[b, q, k]

    return elag, H_mat

end

#for I_MO
function compute_Lagrange2(
    C::Matrix{Float64},
    n::Vector{Float64},
    H::Matrix{Float64},
    I::Array{Float64,4},
    cj12::Matrix{Float64},
    ck12::Matrix{Float64},
    nalpha::Int64,
    nbeta::Int64,
)

    nbf = size(C)[1]
    nbf5 = size(n)[1]
    Cnbf5 = view(C, 1:nbf, 1:nbf5)
    @tullio tmp[m, j] := H[m, nn] * Cnbf5[nn, j]
    @tullio H_mat[i, j] := C[m, i] * tmp[m, j]
    @tullio tmp[m, q, s, l] := Cnbf5[nn, q] * I[m, nn, s, l]
    @tullio tmp2[m, q, r, l] := Cnbf5[s, r] * tmp[m, q, s, l]
    @tullio tmp[m, q, r, t] := Cnbf5[l, t] * tmp2[m, q, r, l]
    @tullio I_MO[p, q, r, t] := C[m, p] * tmp[m, q, r, t]


    elag = zeros(nbf, nbf)

    n_beta = view(n, 1:nbeta)
    n_alpha = view(n, (nalpha+1):nbf5)
    Hmat_nbf5 = view(H_mat, 1:nbf, 1:nbf5)
    grad_nbf5 = view(elag, 1:nbf, 1:nbf5)
    grad_nbeta = view(elag, 1:nbf, 1:nbeta)
    grad_nalpha = view(elag, 1:nbf, (nalpha+1):nbf5)
    I_nb_nb_nb = view(I_MO, 1:nbf, 1:nbeta, 1:nbeta, 1:nbeta)
    I_na_na_na = view(I_MO, 1:nbf, (nalpha+1):nbf5, (nalpha+1):nbf5, (nalpha+1):nbf5)
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

    return elag, H_mat

end


function ENERGY1r(C, n, H, I, cj12, ck12, p)

    if (p.no1 == 0)
        elag, Hmat = compute_Lagrange2(C, n, H, I, cj12, ck12, p.nalpha, p.nbeta)
        E = computeE_elec(Hmat, n, elag, p)
    else
        J, K = computeJKj(C, I, b_mnl, p)

        F = computeF_RC(J, K, n, H, cj12, ck12, p)

        elag = computeLagrange(F, C, p)
    end
    E = computeE_elec(H, C, n, elag, p)

    sumdiff, maxdiff = computeLagrangeConvergency(elag)

    return E, elag, sumdiff, maxdiff


end

function computeF_RC(J, K, n, H, cj12, ck12, p)

    # Matriz de Fock Generalizada
    F = zeros(p.nbf5, p.nbf, p.nbf)

    ini = 0
    if (p.no1 > 1)
        ini = p.no1
    end

    # nH
    @tullio F[i, m, s] += n[i] * H[m, s]

    # nJ
    F_ini_beta = view(F,((ini+1):p.nbeta),:,:)
    n_ini_beta = view(n, (ini+1):p.nbeta)
    J_ini_beta = view(J,((ini+1):p.nbeta),:,:)
    @tullio F_ini_beta[i, m, n] += n_ini_beta[i] * J_ini_beta[i, m, n]
    F_alpha_nbf5 = view(F,((p.nalpha+1):p.nbf5),:,:)
    n_alpha_nbf5 = view(n, (p.nalpha+1):p.nbf5)
    J_alpha_nbf5 = view(J,((p.nalpha+1):p.nbf5),:,:)
    @tullio F_alpha_nbf5[i, m, n] += n_alpha_nbf5[i] * J_alpha_nbf5[i, m, n]

    # C^J J
    cj12_ini_nbf5 = view(cj12, (ini+1):p.nbf5, (ini+1):p.nbf5)
    cj12_ini_nbf5[diagind(cj12_ini_nbf5)] .= 0.0
    @tullio F[i, m, n] += cj12[i, j] * J[j, m, n]

    # -C^K K
    ck12_ini_nbf5 = view(ck12, (ini+1):p.nbf5, (ini+1):p.nbf5)
    ck12_ini_nbf5[diagind(ck12_ini_nbf5)] .= 0.0
    @tullio F[i, m, n] += -ck12[i, j] * K[j, m, n]

    return F

end

function computeLagrange(F, C, p)

    Cnbf5 = view(C, :, 1:p.nbf5)
    Cnoptorb = view(C, :, 1:p.noptorb)

    @tullio G[m, i] := F[i, m, n] * Cnbf5[n, i]

    #Compute Lagrange multipliers
    elag = zeros(p.nbf, p.nbf)
    elag_noptorb_nbf5 = view(elag, 1:p.noptorb, 1:p.nbf5)
    @tullio elag_noptorb_nbf5[i, j] = Cnoptorb[m, i] * G[m, j]

    return elag

end

function computeE_elec(H, C, n, elag, p)
    #EELECTRr
    E = 0
    elag_nbf5 = view(elag, 1:p.nbf5, 1:p.nbf5)
    @tullio E += elag_nbf5[i, i]
    n_beta = view(n, 1:p.nbeta)
    C_beta = view(C, :, 1:p.nbeta)
    @tullio E += n_beta[i] * C_beta[m, i] * H[m, n] * C_beta[n, i]
    if (!p.HighSpin)
        n_beta_alpha = view(n, (p.nbeta+1):p.nalpha)
        C_beta_alpha = view(C, :, (p.nbeta+1):p.nalpha)
        @tullio avx = false E +=
            n_beta_alpha[i] * C_beta_alpha[m, i] * H[m, n] * C_beta_alpha[n, i]
    elseif (p.HighSpin)
        C_beta_alpha = view(C, :, (p.nbeta+1):p.nalpha)
        @tullio avx = false E += 0.5 * C_beta_alpha[m, i] * H[m, n] * C_beta_alpha[n, i]
    end
    n_alpha_nbf5 = view(n, (p.nalpha+1):p.nbf5)
    C_alpha_nbf5 = view(C, :, (p.nalpha+1):p.nbf5)
    @tullio E += n_alpha_nbf5[i] * C_alpha_nbf5[m, i] * H[m, n] * C_alpha_nbf5[n, i]

    return E

end

function computeE_elec(Hmat, n, elag, p)
    #EELECTRr
    E = 0
    elag_nbf5 = view(elag, 1:p.nbf5, 1:p.nbf5)
    @tullio E += elag_nbf5[i, i]
    n_nbf5 = view(n, 1:p.nbf5)
    H_nbf5 = view(Hmat, 1:p.nbf5, 1:p.nbf5)
    @tullio E += n_nbf5[i] * H_nbf5[i, i]

    return E

end

function computeLagrangeConvergency(elag)
    # Convergency
    sumdiff = sum(abs.(elag - Transpose(elag)))
    maxdiff = maximum(abs.(elag - Transpose(elag)))

    return sumdiff, maxdiff

end

function fmiug_scaling(fmiug0, elag, i_ext, nzeros, nbf, noptorb)

    #scaling
    fmiug = zeros(nbf, nbf)
    fmiug_noptorb_noptorb = view(fmiug, 1:noptorb, 1:noptorb)
    if i_ext == 1 && isnothing(fmiug0)
        @tullio fmiug[i, j] = (elag[i, j] + elag[j, i]) / 2
    else
        @tullio fmiug[i, j] = (elag[i, j] - elag[j, i])
        fmiug = tril(fmiug, -1) + Transpose(tril(fmiug, -1))
        for k = 0:(nzeros+9)
            fmiug[10.0^(9-k) .< abs.(fmiug) .< 10.0^(10-k)] .*= 0.1
        end
        fmiug[diagind(fmiug)] .= fmiug0
    end

    return fmiug

end

function fmiug_diis(fk, fmiug, idiis, bdiis, cdiis, maxdiff, p)

    idiis = idiis + 1
    fk[idiis, 1:p.noptorb, 1:p.noptorb] = fmiug[1:p.noptorb, 1:p.noptorb]
    for m = 1:idiis
        bdiis[m, idiis] = 0
        for i = 1:p.noptorb
            for j = 1:i
                bdiis[m, idiis] += fk[m, i, j] * fk[idiis, j, i]
            end
        end
        bdiis[idiis, m] = bdiis[m, idiis]
        bdiis[m, idiis+1] = -1
        bdiis[idiis+1, m] = -1
    end
    bdiis[idiis+1, idiis+1] = 0

    if idiis >= p.ndiis
        cdiis = zeros(idiis + 1)
        cdiis[1:idiis] .= 0
        cdiis[idiis+1] .= -1
        x = bdiis[1:(idiis+1), 1:(idiis+1)] \ cdiis[1:(idiis+1)]

        for i = 1:p.noptorb
            for j = 1:i
                fmiug[i, j] = 0
                for k = 1:(idiis+1)
                    fmiug[i, j] += x[k] * fk[k, i, j]
                end
                fmiug[j, i] = fmiug[i, j]
            end
        end
    end

    if p.perdiis
        idiis = 0
    end

    return fk, fmiug, idiis, bdiis
end

function compute_E_nuc(bset, p)

    E_nuc = 0.0
    for i = 1:p.natoms
        for j = (i+1):p.natoms
            E_nuc +=
                bset.atoms[i].Z * bset.atoms[j].Z /
                (norm(bset.atoms[i].xyz - bset.atoms[j].xyz) * 1.88973)
        end
    end

    return E_nuc

end

function Z_to_symbol(Z)
    dict = Dict(
        1 => "H",
        2 => "He",
        3 => "Li",
        4 => "Be",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        10 => "Ne",
        11 => "Na",
        12 => "Mg",
        13 => "Al",
        14 => "Si",
        15 => "P",
        16 => "S",
        17 => "Cl",
        18 => "Ar",
        19 => "K",
        20 => "Ca",
        21 => "Sc",
        22 => "Ti",
        23 => "V",
        24 => "Cr",
        25 => "Mn",
        26 => "Fe",
        27 => "Co",
        28 => "Ni",
        29 => "Cu",
        30 => "Zn",
        31 => "Ga",
        32 => "Ge",
        33 => "As",
        34 => "Se",
        35 => "Br",
        36 => "Kr",
        37 => "Rb",
        38 => "Sr",
        39 => "Y",
        40 => "Zr",
        41 => "Nb",
        42 => "Mo",
        43 => "Tc",
        44 => "Ru",
        45 => "Rh",
        46 => "Pd",
        47 => "Ag",
        48 => "Cd",
        49 => "In",
        50 => "Sn",
        51 => "Sb",
        52 => "Te",
        53 => "I",
        54 => "Xe",
        55 => "Cs",
        56 => "Ba",
        57 => "La",
        58 => "Ce",
        59 => "Pr",
        60 => "Nd",
        61 => "Pm",
        62 => "Sm",
        63 => "Eu",
        64 => "Gd",
        65 => "Tb",
        66 => "Dy",
        67 => "Ho",
        68 => "Er",
        69 => "Tm",
        70 => "Yb",
        71 => "Lu",
        72 => "Hf",
        73 => "Ta",
        74 => "W",
        75 => "Re",
        76 => "Os",
        77 => "Ir",
        78 => "Pt",
        79 => "Au",
        80 => "Hg",
        81 => "Tl",
        82 => "Pb",
        83 => "Bi",
        84 => "Po",
        85 => "At",
        86 => "Rn",
    )

    return dict[Z]
end

function check_ortho(C, S, p)

    # Revisa ortonormalidad
    orthonormality = true
    CTSC = C' * S * C
    ortho_deviation = abs.(CTSC - I)
    if (maximum(ortho_deviation) > 10^-6)
        orthonormality = false
    end
    if !orthonormality
        @printf(
            "Orthonormality violations %i, Maximum Violation %f\n",
            sum(ortho_deviation .> 10^-6),
            maximum(ortho_deviation)
        )
        println("Trying to orthonormalize")
        C = orthonormalize(C, S, p)
        C = check_ortho(C, S, p)
    else
        println("No violations of the orthonormality")
    end
    for j = 1:p.nbf
        #Obtiene el índice del coeficiente con mayor valor absoluto del MO
        idxmaxabsval = 1
        for i = 1:p.nbf
            if (abs(C[i, j]) > abs(C[idxmaxabsval, j]))
                idxmaxabsval = i
            end
        end
        # Ajusta el signo del MO
        C[1:p.nbf, j] = sign(C[idxmaxabsval, j]) * C[1:p.nbf, j]
    end

    return C

end

function orthonormalize(C, S, p)

    evals, evecs = eigen(S)
    for i = 1:size(evals)[1]
        if (evals[i] < 0.0)
            evals[i] = 0.0
        else
            evals[i] = 1 / sqrt(evals[i])
        end
    end

    @tullio S_12[i, j] := evecs[i, j] * evals[j]

    @tullio Cnew[i, j] := S[i, k] * C[k, j]

    @tullio Cnew2[i, j] := S_12[k, i] * Cnew[k, j]

    for i = 1:size(Cnew2)[2]
        Cnew2[1:p.nbf, i] = Cnew2[1:p.nbf, i] / norm(Cnew2[1:p.nbf, i])
        for j = (i+1):size(Cnew2)[1]
            val = -sum(Cnew2[1:p.nbf, i] .* Cnew2[1:p.nbf, j])
            Cnew2[1:p.nbf, j] = Cnew2[1:p.nbf, j] + val * Cnew2[1:p.nbf, i]
        end
    end

    @tullio C[i, j] = S_12[i, k] * Cnew2[k, j]

    return C

end

function rotate_orbital(y, C, p)

    ynew = zeros(p.nbf, p.nbf)

    n = 1
    for i = 1:p.nbf5
        for j = (i+1):p.nbf
            ynew[i, j] = y[n]
            ynew[j, i] = -y[n]
            n += 1
        end
    end

    #tmp = [i<j&&i<=p.nbf5 ? y[Int64((2*p.nbf*i - i^2 - i)/2 + j - p.nbf)] : 0 for i in 1:p.nbf, j in 1:p.nbf]
    #ynew = tmp .- tmp'

    U = exp(ynew)
    U = real.(U)
    @tullio Cnew[m, p] := C[m, r] * U[r, p]

    return Cnew

end

function n_to_gammas_ebi(n)
    """Transform n to gammas in the ebi encoding

    x_p = erf^-1 (2n_p - 1)

    """

    nv = size(n)[1]
    gamma = zeros(nv)
    for i = 1:nv
        gamma[i] = erfinv(2 * n[i] - 1)
    end
    return gamma

end

function n_to_gammas_softmax(n, p)

    gamma = zeros(p.nv)

    for i = 1:p.ndoc

        ll = p.no1 + p.ndns + (p.ndoc - i) + 1
        ul = ll + p.ndoc * (p.ncwo - 1)
        n_pi = @view n[ll:p.ndoc:ul]

        llg = (p.ndoc - i) + 1
        ulg = llg + p.ndoc * (p.ncwo - 1)
        gamma_pi = @view gamma[llg:p.ndoc:ulg]

        for j = 1:p.ncwo
            gamma_pi[j] = log(max(n_pi[j], 1e-16)) - log(n[i])
        end

    end

    return gamma

end

function n_to_gammas_trigonometric(n, p)
    gamma = zeros(p.nv)
    for i = 1:p.ndoc
        idx = p.no1 + i
        gamma[i] = acos(sqrt(2.0 * n[idx] - 1.0))
        prefactor = max(1 - n[idx], 1e-14)

        ll_n = p.no1 + p.ndns + (p.ndoc - i) + 1
        ul_n = ll_n + p.ndoc * (p.ncwo - 1)
        n_pi = @view n[ll_n:p.ndoc:ul_n]

        ll_gamma = p.ndoc + (p.ndoc - i) + 1
        ul_gamma = ll_gamma + p.ndoc * (p.ncwo - 2)
        gamma_pi = @view gamma[ll_gamma:p.ndoc:ul_gamma]

        for j = 1:(p.ncwo-1)
            gamma_pi[j] = asin(sqrt(n_pi[j] / prefactor))
            prefactor = prefactor * (cos(gamma_pi[j]))^2
        end
    end
    return gamma
end

function order_occupations_softmax(old_C, old_gamma, p)

    C = copy(old_C)
    gamma = copy(old_gamma)

    if any(gamma .> 0)
        #Sort ndoc subspaces
        C_tmp = zeros(p.nbf, 1 + p.ncwo)
        for i = 1:p.ndoc
            old_ll = p.no1 + p.ndns + (p.ndoc - i) + 1
            old_ul = old_ll + p.ndoc * (p.ncwo - 1)

            C_tmp[:, 1] = old_C[:, p.no1+i]
            C_tmp[:, 2:end] = old_C[:, old_ll:p.ndoc:old_ul]

            old_ll_x = (p.ndoc - i) + 1
            old_ul_x = old_ll_x + p.ndoc * (p.ncwo - 1)

            gamma_tmp = vcat([0], gamma[old_ll_x:p.ndoc:old_ul_x])
            sort_idx = sortperm(gamma_tmp)[end:-1:1]

            gamma_tmp = gamma_tmp[sort_idx] .- maximum(gamma_tmp)
            C_tmp = C_tmp[:, sort_idx]

            gamma[old_ll_x:p.ndoc:old_ul_x] = gamma_tmp[2:end]
            C[:, p.no1+i] = C_tmp[:, 1]
            C[:, old_ll:p.ndoc:old_ul] = C_tmp[:, 2:end]
        end
    end
    return C, gamma
end
