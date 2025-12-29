function PNOFi_selector(n, p)
    if p.ipnof == 4
        cj12, ck12 =
            CJCKD4(n, p.no1, p.ndoc, p.nsoc, p.nbeta, p.nalpha, p.ndns, p.ncwo, p.MSpin)
    elseif p.ipnof == 5
        cj12, ck12 =
            CJCKD5(n, p.no1, p.ndoc, p.nsoc, p.nbeta, p.nalpha, p.ndns, p.ncwo, p.MSpin)
    elseif p.ipnof == 7
        cj12, ck12 = CJCKD7(
            n,
            p.ista,
            p.no1,
            p.ndoc,
            p.nsoc,
            p.nbeta,
            p.nalpha,
            p.ndns,
            p.ncwo,
            p.MSpin,
        )
    elseif p.ipnof == 8
        cj12, ck12 = CJCKD8(
            n,
            p.no1,
            p.ndoc,
            p.nsoc,
            p.nbeta,
            p.nalpha,
            p.ndns,
            p.ncwo,
            p.MSpin,
            p.h_cut,
            p.ista,
        )
    elseif p.ipnof == 9
        cj12, ck12 = CJCKD9(
            n,
            p.no1,
            p.ndoc,
            p.nsoc,
            p.nbeta,
            p.nalpha,
            p.ndns,
            p.ncwo,
            p.MSpin,
            p.h_cut,
        )
    end
    cj12[diagind(cj12)] .= 0
    ck12[diagind(ck12)] .= 0

    return cj12, ck12
end

function der_PNOFi_selector(n, dn_dgamma, p)
    if p.ipnof == 4
        Dcj12r, Dck12r =
            der_CJCKD4(n, dn_dgamma, p.no1, p.ndoc, p.nalpha, p.nv, p.nbf5, p.ndns, p.ncwo)
    elseif p.ipnof == 5
        Dcj12r, Dck12r =
            der_CJCKD5(n, dn_dgamma, p.no1, p.ndoc, p.nalpha, p.nv, p.nbf5, p.ndns, p.ncwo)
    elseif p.ipnof == 7
        Dcj12r, Dck12r = der_CJCKD7(
            n,
            p.ista,
            dn_dgamma,
            p.no1,
            p.ndoc,
            p.nalpha,
            p.nv,
            p.nbf5,
            p.ndns,
            p.ncwo,
        )
    elseif p.ipnof == 8
        Dcj12r, Dck12r = der_CJCKD8(
            n,
            dn_dgamma,
            p.no1,
            p.ndoc,
            p.nalpha,
            p.nbeta,
            p.nv,
            p.nbf5,
            p.ndns,
            p.ncwo,
            p.MSpin,
            p.nsoc,
            p.h_cut,
            p.ista,
        )
    elseif p.ipnof == 9
        Dcj12r, Dck12r = der_CJCKD9(
            n,
            dn_dgamma,
            p.no1,
            p.ndoc,
            p.nalpha,
            p.nbeta,
            p.nv,
            p.nbf5,
            p.ndns,
            p.ncwo,
            p.MSpin,
            p.nsoc,
            p.h_cut,
        )
    end
    for i = 1:p.nv
        for j = 1:p.nbf5
            Dcj12r[j, j, i] = 0
            Dck12r[j, j, i] = 0
        end
    end

    return Dcj12r, Dck12r
end

function CJCKD4(n, no1, ndoc, nsoc, nbeta, nalpha, ndns, ncwo, MSpin)
    """PNOF4 coefficients C^J and C^K that multiply J and K integrals

    E = 2sum_p n_pH_p + sum_{pq}' C^J_{pq}J_{qp} - sum_{pq}' C^K_{pq}K_{qp}

    T = (1-S_F)/S_F
    R_pq = h_p n_q/S_F; p in occ, q in vir

    Delta_{pq} = begin{cases}
                 h_p h_q   & if p in occ, q in occ
                 T h_p n_q & if p in occ, q in vir
                 T n_p h_q & if p in vir, q in occ
                 n_p n_q   & if p in vir, q in vir
               end{cases}

    Pi_{pq} = begin{cases}
                -sqrt{h_p h_q}                     & if p in occ, q in occ
                -sqrt{R_pq} sqrt{n_p - n_q + R_pq} & if p in occ, q in vir
                -sqrt{R_qp} sqrt{n_q - n_p + R_qp} & if p in vir, q in occ
                 sqrt{n_p n_q}                     & if p in vir, q in vir
               end{cases}

    C^J = 2(n_pn_q - Delta_pq)
    C^K = n_pn_q - Delta_pq - Pi_pq
    """

    nbf5 = size(n)[1]

    h = 1 .- n
    S_F = sum(h[1:ndoc])

    T = (1 - S_F) / S_F
    R = (h[1:ndoc] * n[(ndoc+1):end]') / S_F

    Delta = zeros(nbf5, nbf5)
    Delta[1:ndoc, 1:ndoc] = h[1:ndoc] * h[1:ndoc]'
    Delta[1:ndoc, (ndoc+1):end] = (1 - S_F) / S_F * h[1:ndoc] * n[(ndoc+1):end]'
    Delta[(ndoc+1):end, 1:ndoc] = (1 - S_F) / S_F * n[(ndoc+1):end] * h[1:ndoc]'
    Delta[(ndoc+1):end, (ndoc+1):end] = n[(ndoc+1):end] * n[(ndoc+1):end]'

    Pi = zeros(nbf5, nbf5)
    Pi[1:ndoc, 1:ndoc] = -sqrt.(h[1:ndoc] * h[1:ndoc]')
    Pi[1:ndoc, (ndoc+1):end] =
        -sqrt.((h[1:ndoc] * n[(ndoc+1):end]') / S_F) .*
        sqrt.((n[1:ndoc] .- n[(ndoc+1):end]') .+ R)
    Pi[(ndoc+1):end, 1:ndoc] =
        -sqrt.((n[(ndoc+1):end] * h[1:ndoc]') / S_F) .*
        sqrt.(-(n[(ndoc+1):end] .- n[1:ndoc]') .+ R')
    Pi[(ndoc+1):end, (ndoc+1):end] = sqrt.(n[(ndoc+1):end] * n[(ndoc+1):end]')

    cj12 = 2 * ((n * n') - Delta)
    ck12 = (n * n') - Delta - Pi

    return cj12, ck12
end

function der_CJCKD4(n, dn_dgamma, no1, ndoc, nalpha, nv, nbf5, ndns, ncwo)
    """PNOF4 coefficients DC^J and DC^K that multiply J and K integrals

    dE/dgamma_g = 2sum_p dn_p/dgamma_g H_p 
                  +sum_{pq} dC^J_{pq}/dgamma_g J_{qp}
                  -sum_{pq} dC^K_{pq}/dgamma_g K_{qp}

    T = (1-S_F)/S_F
    R_pq = h_p n_q/S_F; p in occ, q in vir

    dS_F/dx_g = - sum_p^F dn_p/dx_g
    dT/dx_g = -(dS_F/dx_g)/S_F^2
    dR_pq/dx_g = [(dh_p/dx_g n_q + h_p dn_q/dx_g) S_F - (h_p n_q) dS_F/dx_g ] / S_F^2

    dDelta_{pq}/dx_g = begin{cases}
                 dh_p/dx_g h_q + h_p dh_q/dx_g                       & if p in occ, q in occ
                 dT/dx_g h_p n_q + T dh_p/dx_g n_q + T h_p dn_q/dx_g & if p in occ, q in vir
                 dT/dx_g n_p h_q + T dn_p/dx_g h_q + T n_p dh_q/dx_g & if p in vir, q in occ
                 dn_p/dx_g n_q + n_p dn_q/dx_g                       & if p in vir, q in vir
               end{cases}

    dPi_{pq}/dx_g = begin{cases}
      -1/2*(1/sqrt{h_p} dh_p/dx_g sqrt{h_q} + sqrt{h_p} 1/sqrt{h_q} dh_q/dx_g )  & if p in occ, q in occ
      -((1/2sqrt{R_pq}) dR_pq/dx_g sqrt{n_p - n_q + R_pq} + sqrt{R_pq}(dn_p/dx_g - dn_q/dx_g + dR_pq/dx_g)/(2sqrt{n_q-n_p+R_pq})) & if p in occ, q in vir
      -((1/2sqrt{R_qp}) dR_qp/dx_g sqrt{n_q - n_p + R_pq} + sqrt{R_qp}(dn_q/dx_g - dn_p/dx_g + dR_qp/dx_g)/(2sqrt{n_p-n_q+R_qp})) & if p in vir, q in occ
      1/2*(1/sqrt{n_p} dn_p/dx_g sqrt{n_q} + sqrt{n_p} 1/sqrt{n_q} dn_q/dx_g )  & if p in vir, q in vir
      end{cases}

    """

    nbf5 = size(n)[1]

    h = 1 .- n
    dh_dgamma = -dn_dgamma

    S_F = max(sum(h[1:ndoc]), 1e-7)
    D_S_F = zeros(nv)
    for k = 1:nv
        D_S_F[k] = sum(dh_dgamma[1:ndoc, k])
    end
    T = (1 - S_F) / S_F
    D_T = zeros(nv)
    for k = 1:nv
        D_T[k] = (-D_S_F[k] / S_F^2)
    end

    R = (h[1:ndoc] * n[(ndoc+1):end]') / S_F
    D_R = zeros(ndoc, nbf5 - ndoc, nv)
    for k = 1:nv
        D_R[1:end, 1:end, k] =
            (
                (
                    (dh_dgamma[1:ndoc, k] * n[(ndoc+1):end]') .+
                    (h[1:ndoc] * dn_dgamma[(ndoc+1):end, k]')
                ) * S_F .- (h[1:ndoc] * n[(ndoc+1):end]') * D_S_F[k]
            ) / S_F^2
    end

    D_Delta = zeros(nbf5, nbf5, nv)
    for k = 1:nv
        D_Delta[1:ndoc, 1:ndoc, k] = dh_dgamma[1:ndoc, k] * h[1:ndoc]'
        D_Delta[1:ndoc, (ndoc+1):end, k] =
            D_T[k] .* (h[1:ndoc] * n[(ndoc+1):end]') .+
            T .* (dh_dgamma[1:ndoc, k] * n[(ndoc+1):end]')
        D_Delta[(ndoc+1):end, 1:ndoc, k] = T * (dn_dgamma[(ndoc+1):end, k] * h[1:ndoc]')
        D_Delta[(ndoc+1):end, (ndoc+1):end, k] =
            dn_dgamma[(ndoc+1):end, k] * n[(ndoc+1):end]'
    end

    D_Pi = zeros(nbf5, nbf5, nv)
    for k = 1:nv
        D_Pi[1:ndoc, 1:ndoc, k] =
            -(1 ./ max.(2 * sqrt.(h[1:ndoc]), 1e-7) .* dh_dgamma[1:ndoc, k]) *
            sqrt.(h[1:ndoc])'
        D_Pi[1:ndoc, (ndoc+1):end, k] = -(
            1 ./ max.(2 * sqrt.(R), 1e-7) .* D_R[1:end, 1:end, k] .*
            sqrt.((n[1:ndoc] .- n[(ndoc+1):end]') .+ R) .+
            sqrt.(R) .* (
                (dn_dgamma[1:ndoc, k] .- dn_dgamma[(ndoc+1):end, k]') .+
                D_R[1:end, 1:end, k]
            ) ./ max.(2 * sqrt.((n[1:ndoc] * n[(ndoc+1):end]') .+ R), 1e-7)
        )
        D_Pi[(ndoc+1):end, (ndoc+1):end, k] =
            (1 ./ max.(2 * sqrt.(n[(ndoc+1):end]), 1e-7) .* dn_dgamma[(ndoc+1):end, k]) *
            (sqrt.(n[(ndoc+1):end]))'
    end

    Dcj12r = zeros(nbf5, nbf5, nv)
    Dck12r = zeros(nbf5, nbf5, nv)
    for k = 1:nv
        Dcj12r[1:end, 1:end, k] =
            2 * ((dn_dgamma[1:end, k] * n') - D_Delta[1:end, 1:end, k])
        Dck12r[1:end, 1:end, k] =
            (dn_dgamma[1:end, k] * n') - D_Delta[1:end, 1:end, k] - D_Pi[1:end, 1:end, k]
    end

    return Dcj12r, Dck12r
end

function CJCKD5(n, no1, ndoc, nsoc, nbeta, nalpha, ndns, ncwo, MSpin)

    # Interpair Electron correlation #

    @tullio cj12[i, j] := 2 * n[i] * n[j]
    @tullio ck12[i, j] := n[i] * n[j]

    # Intrapair Electron Correlation #

    if MSpin == 0 && nsoc > 1
        n_beta_alpha = view(n, (nbeta+1):nalpha)
        ck12_beta_alpha = view(ck12, (nbeta+1):nalpha, (nbeta+1):nalpha)
        @tullio ck12_beta_alpha[i, j] = 2 * n_beta_alpha[i] * n_beta_alpha[j]
    end

    for l = 1:ndoc
        ldx = no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - l) + 1
        ul = ll + ndoc * (ncwo - 1)

        cj12[ldx, ll:ndoc:ul] .= 0
        cj12[ll:ndoc:ul, ldx] .= 0

        cj12[ll:ndoc:ul, ll:ndoc:ul] .= 0

        ck12[ldx, ll:ndoc:ul] .= sqrt.(n[ldx] * n[ll:ndoc:ul])
        ck12[ll:ndoc:ul, ldx] .= sqrt.(n[ldx] * n[ll:ndoc:ul])

        ck12_ww = view(ck12, ll:ndoc:ul, ll:ndoc:ul)
        n_ww = view(n, ll:ndoc:ul)
        @tullio ck12_ww[i, j] = -sqrt(n_ww[i] * n_ww[j])
    end

    return cj12, ck12

end

function der_CJCKD5(n, dn_dgamma, no1, ndoc, nalpha, nv, nbf5, ndns, ncwo)

    # Interpair Electron correlation #

    @tullio Dcj12r[i, j, k] := 2 * dn_dgamma[i, k] * n[j]
    @tullio Dck12r[i, j, k] := dn_dgamma[i, k] * n[j]

    # Intrapair Electron Correlation

    for l = 1:ndoc
        ldx = no1 + l

        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - l) + 1
        ul = ll + ndoc * (ncwo - 1)

        Dcj12r[ldx, ll:ndoc:ul, 1:nv] .= 0
        Dcj12r[ll:ndoc:ul, ldx, 1:nv] .= 0

        Dcj12r[ll:ndoc:ul, ll:ndoc:ul, 1:nv] .= 0

        a = max(n[ldx], 10^-15)
        b = n[ll:ndoc:ul]
        b[b .< 10^-15] .= 10^-15

        Dck12r_occ_cwo = view(Dck12r, ldx, ll:ndoc:ul, 1:nv)
        Dck12r_cwo_occ = view(Dck12r, ll:ndoc:ul, ldx, 1:nv)
        Dck12r_cwo_cwo = view(Dck12r, ll:ndoc:ul, ll:ndoc:ul, 1:nv)
        n_cwo = view(n, ll:ndoc:ul)
        n_occ = n[ldx]
        dn_dgamma_occ = view(dn_dgamma, ldx, 1:nv)
        dn_dgamma_cwo = view(dn_dgamma, ll:ndoc:ul, 1:nv)
        @tullio Dck12r_occ_cwo[i, j] =
            1 / 2 * 1 / sqrt(a) * dn_dgamma_occ[j] * sqrt(n_cwo[i])
        @tullio Dck12r_cwo_occ[i, j] =
            1 / 2 * 1 / sqrt(b[i]) * dn_dgamma_cwo[i, j] * sqrt(n_occ)
        @tullio Dck12r_cwo_cwo[i, j, k] =
            -1 / 2 * 1 / sqrt(b[i]) * dn_dgamma_cwo[i, k] * sqrt(n_cwo[j])

    end
    return Dcj12r, Dck12r

end

function CJCKD7(
    n::Vector{Float64},
    ista::Int64,
    no1::Int64,
    ndoc::Int64,
    nsoc::Int64,
    nbeta::Int64,
    nalpha::Int64,
    ndns::Int64,
    ncwo::Int64,
    MSpin::Int64,
)

    if ista == 0
        fi = n .* (1 .- n)
        fi[fi .<= 0] .= 0
        fi = sqrt.(fi)
    else
        fi = 2 * n .* (1 .- n)
    end

    # Interpair Electron correlation #

    @tullio grad=false cj12[i, j] := 2 * n[i] * n[j]
    @tullio grad=false ck12[i, j] := n[i] * n[j] + fi[i] * fi[j]

    # Intrapair Electron Correlation

    if nsoc > 1
        n_beta_alpha = view(n, (nbeta+1):nalpha)
        ck12_beta_alpha = view(ck12, (nbeta+1):nalpha, (nbeta+1):nalpha)
        @tullio grad=false ck12_beta_alpha[i, j] = 2 * n_beta_alpha[i] * n_beta_alpha[j]
    end

    for l = 1:ndoc
        ldx = no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - l) + 1
        ul = ll + ndoc * (ncwo - 1)

        cj12[ldx, ll:ndoc:ul] .= 0
        cj12[ll:ndoc:ul, ldx] .= 0

        cj12[ll:ndoc:ul, ll:ndoc:ul] .= 0

        ck12[ldx, ll:ndoc:ul] .= sqrt.(n[ldx] * n[ll:ndoc:ul])
        ck12[ll:ndoc:ul, ldx] .= sqrt.(n[ldx] * n[ll:ndoc:ul])

        ck12_ww = view(ck12, ll:ndoc:ul, ll:ndoc:ul)
        n_ww = view(n, ll:ndoc:ul)
        @tullio grad=false ck12_ww[i, j] = -sqrt(n_ww[i] * n_ww[j])
    end

    return cj12, ck12

end

function der_CJCKD7(
    n::Vector{Float64},
    ista::Int64,
    dn_dgamma::Matrix{Float64},
    no1::Int64,
    ndoc::Int64,
    nalpha::Int64,
    nv::Int64,
    nbf5::Int64,
    ndns::Int64,
    ncwo::Int64,
)

    if ista == 0
        fi = n .* (1 .- n)
        fi[fi .<= 0] .= 0
        fi = sqrt.(fi)
    else
        fi = 2 * n .* (1 .- n)
    end

    dfi_dgamma = zeros(nbf5, nv)
    for i = (no1+1):nbf5
        a = max(fi[i], 10^-15)
        for k = 1:nv
            if ista == 0
                dfi_dgamma[i, k] = 1 / (2 * a) * (1 - 2 * n[i]) * dn_dgamma[i, k]
            else
                dfi_dgamma[i, k] = 2 * (1 - 2 * n[i]) * dn_dgamma[i, k]
            end
        end
    end

    # Interpair Electron correlation #

    @tullio grad=false Dcj12r[i, j, k] := 2 * dn_dgamma[i, k] * n[j]
    @tullio grad=false Dck12r[i, j, k] := dn_dgamma[i, k] * n[j] + dfi_dgamma[i, k] * fi[j]

    # Intrapair Electron Correlation

    for l = 1:ndoc
        ldx = no1 + l

        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - l) + 1
        ul = ll + ndoc * (ncwo - 1)

        Dcj12r[ldx, ll:ndoc:ul, 1:nv] .= 0
        Dcj12r[ll:ndoc:ul, ldx, 1:nv] .= 0

        Dcj12r[ll:ndoc:ul, ll:ndoc:ul, 1:nv] .= 0

        a = max(n[ldx], 10^-15)
        b = n[ll:ndoc:ul]
        b[b .< 10^-15] .= 10^-15

        Dck12r_occ_cwo = view(Dck12r, ldx, ll:ndoc:ul, 1:nv)
        Dck12r_cwo_occ = view(Dck12r, ll:ndoc:ul, ldx, 1:nv)
        Dck12r_cwo_cwo = view(Dck12r, ll:ndoc:ul, ll:ndoc:ul, 1:nv)
        n_cwo = view(n, ll:ndoc:ul)
        n_occ = n[ldx]
        dn_dgamma_occ = view(dn_dgamma, ldx, 1:nv)
        dn_dgamma_cwo = view(dn_dgamma, ll:ndoc:ul, 1:nv)
        @tullio grad=false Dck12r_occ_cwo[i, j] =
            1 / 2 * 1 / sqrt(a) * dn_dgamma_occ[j] * sqrt(n_cwo[i])
        @tullio grad=false Dck12r_cwo_occ[i, j] =
            1 / 2 * 1 / sqrt(b[i]) * dn_dgamma_cwo[i, j] * sqrt(n_occ)
        @tullio grad=false Dck12r_cwo_cwo[i, j, k] =
            -1 / 2 * 1 / sqrt(b[i]) * dn_dgamma_cwo[i, k] * sqrt(n_cwo[j])

    end
    return Dcj12r, Dck12r

end

function CJCKD8(n, no1, ndoc, nsoc, nbeta, nalpha, ndns, ncwo, MSpin, h_cut, ista)

    # ista = 0 GNOF
    # ista = 1 GNOFm (with N<F static interpair interactions)
    # ista = 2 GNOFs (with N<F static interpair interactions and phi = 0.9 sqrt(n*(1-n)) )
    # ista = 2 GNOFmsta (with N<F static interpair interactions and phi = a (n*(1-n))^0.55 )
    # ista = 3 GNOFmodulated (with N<F static interpair interactions and phi = exp(-(n-0.5)^2) sqrt(n(1-n)) )

    nbf5 = size(n)[1]

    #h_cut = p.h_cut#0.02*sqrt(2.0)
    n_d = zeros(nbf5)

    c = 1.0
    if ista == 5
        c = 0.8
    end

    for i = 1:ndoc
        idx = no1 + i
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - i) + 1
        ul = ll + ndoc * (ncwo - 1)

        h = 1.0 - n[idx]
        coc = h / h_cut
        arg = -coc^2
        F = c*exp(arg)  # ! Hd/Hole
        n_d[idx] = n[idx] * F
        n_d[ll:ndoc:ul] = n[ll:ndoc:ul] * F  # ROd = RO*Hd/Hole

    end

    n_d12 = sqrt.(n_d)
    fi = n .* (1 .- n)
    fi[fi .<= 0] .= 0

    if ista == 0 || ista == 1 || ista == 5
        fi = sqrt.(fi)
    elseif ista == 2
        alpha = 0.9
        fi = alpha * sqrt.(fi)
    elseif ista == 3
        alpha = 0.55
        aa = 2.0^(2.0*alpha-1)
        aa = aa ^ (1/alpha)
        fi = aa * fi
        fi = fi .^ alpha
    elseif ista == 4
        fi = sqrt.(fi) .* exp.(- (n .- 0.5) .^ 2)
    end

    # Interpair Electron correlation #

    @tullio cj12[i, j] := 2 * n[i] * n[j]
    if ista == 0
        @tullio ck12[i, j] := n[i] * n[j]

        ck12_beta_cwo = view(ck12, (no1+1):nbeta, (nalpha+1):nbf5)
        ck12_cwo_beta = view(ck12, (nalpha+1):nbf5, (no1+1):nbeta)
        ck12_cwo_cwo = view(ck12, (nalpha+1):nbf5, (nalpha+1):nbf5)
        fi_beta = view(fi, (no1+1):nbeta)
        fi_cwo = view(fi, (nalpha+1):nbf5)
        @tullio ck12_beta_cwo[i, j] += fi_beta[i] * fi_cwo[j]
        @tullio ck12_cwo_beta[i, j] += fi_cwo[i] * fi_beta[j]
        @tullio ck12_cwo_cwo[i, j] += fi_cwo[i] * fi_cwo[j]
        if nsoc > 0
            ck12[(no1+1):nbeta, (nbeta+1):nalpha] .+= 0.5 * fi[(no1+1):nbeta] * 0.5
            ck12[(nbeta+1):nalpha, (no1+1):nbeta] .+= 0.5 * 0.5 * fi[(no1+1):nbeta]'
            ck12[(nbeta+1):nalpha, (nalpha+1):nbf5] .+= 0.5 * fi[(nalpha+1):nbf5]'
            ck12[(nalpha+1):nbf5, (nbeta+1):nalpha] .+= fi[(nalpha+1):nbf5] * 0.5
        end
        if nsoc > 1
            ck12[(nbeta+1):nalpha, (nbeta+1):nalpha] .= 0.5
        end
    elseif ista > 0
        @tullio ck12[i, j] := n[i] * n[j] + fi[i] * fi[j]
        ck12_beta_cwo = view(ck12, (no1+1):nbeta, (nalpha+1):nbf5)
        ck12_cwo_beta = view(ck12, (nalpha+1):nbf5, (no1+1):nbeta)
        ck12_cwo_cwo = view(ck12, (nalpha+1):nbf5, (nalpha+1):nbf5)
    end

    # Intrapair Electron Correlation

    n_d12_beta = view(n_d12, (no1+1):nbeta)
    n_d12_cwo = view(n_d12, (nalpha+1):nbf5)
    n_d_beta = view(n_d, (no1+1):nbeta)
    n_d_cwo = view(n_d, (nalpha+1):nbf5)
    @tullio ck12_beta_cwo[i, j] += n_d12_beta[i] * n_d12_cwo[j] - n_d_beta[i] * n_d_cwo[j]
    @tullio ck12_cwo_beta[i, j] += n_d12_cwo[i] * n_d12_beta[j] - n_d_cwo[i] * n_d_beta[j]
    @tullio ck12_cwo_cwo[i, j] += -n_d12_cwo[i] * n_d12_cwo[j] - n_d_cwo[i] * n_d_cwo[j]

    for l = 1:ndoc
        ldx = no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - l) + 1
        ul = ll + ndoc * (ncwo - 1)

        cj12[ldx, ll:ndoc:ul] .= 0
        cj12[ll:ndoc:ul, ldx] .= 0

        cj12[ll:ndoc:ul, ll:ndoc:ul] .= 0

        ck12[ldx, ll:ndoc:ul] .= sqrt.(n[ldx] * n[ll:ndoc:ul])
        ck12[ll:ndoc:ul, ldx] .= sqrt.(n[ldx] * n[ll:ndoc:ul])

        ck12_ww = view(ck12, ll:ndoc:ul, ll:ndoc:ul)
        n_ww = view(n, ll:ndoc:ul)
        @tullio ck12_ww[i, j] = -sqrt(n_ww[i] * n_ww[j])
    end

    return cj12, ck12

end

function der_CJCKD8(
    n,
    dn_dgamma,
    no1,
    ndoc,
    nalpha,
    nbeta,
    nv,
    nbf5,
    ndns,
    ncwo,
    MSpin,
    nsoc,
    h_cut,
    ista,
)

    nbf5 = size(n)[1]

    #h_cut = p.h_cut#0.02*sqrt(2.0)
    n_d = zeros(nbf5)
    dn_d_dgamma = zeros(nbf5, nv)
    dn_d12_dgamma = zeros(nbf5, nv)

    c = 1.0
    if ista == 5
        c = 0.8
    end
    for i = 1:ndoc
        idx = no1 + i
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - i) + 1
        ul = ll + ndoc * (ncwo - 1)

        h_idx = 1.0 - n[idx]
        coc = h_idx / h_cut
        arg = -coc^2
        F_idx = c * exp(arg)                # Hd/Hole
        n_d[idx] = n[idx] * F_idx
        n_d[ll:ndoc:ul] = n[ll:ndoc:ul] * F_idx      # n_d = RO*Hd/Hole
        dn_d_dgamma[idx, 1:nv] .=
            F_idx * dn_dgamma[idx, 1:nv] .* (1 - n[idx] * (-2 * coc / h_cut))
        dn_d_dgamma[ll:ndoc:ul, 1:nv] .= F_idx * dn_dgamma[ll:ndoc:ul, 1:nv]
        dn_d_dgamma_cwo_v = view(dn_d_dgamma, ll:ndoc:ul, 1:nv)
        dn_dgamma_v = view(dn_dgamma, idx, 1:nv)
        n_cwo = view(n, ll:ndoc:ul)
        @tullio dn_d_dgamma_cwo_v[i, j] +=
            F_idx * (2 * coc / h_cut) * n_cwo[i] * dn_dgamma_v[j]
    end

    n_d12 = sqrt.(n_d)
    dn_d12_dgamma = 0.5 * dn_d_dgamma ./ max.(n_d12, 10^-15)

    fi = n .* (1 .- n)
    fi[fi .<= 0] .= 0

    dfi_dgamma = zeros(nbf5, nv)
    if ista == 0 || ista == 1 || ista == 5
        fi = sqrt.(fi)
        for i = (no1+1):nbf5
            a = max(fi[i], 10^-15)
            for k = 1:nv
                dfi_dgamma[i, k] = 1 / (2 * a) * (1 - 2 * n[i]) * dn_dgamma[i, k]
            end
        end
    elseif ista == 2
        alpha = 0.9
        fi = alpha * sqrt.(fi)
        for i = (no1+1):nbf5
            a = max(fi[i], 10^-15)
            for k = 1:nv
                dfi_dgamma[i, k] = alpha^2 / (2 * a) * (1 - 2 * n[i]) * dn_dgamma[i, k]
            end
        end
    elseif ista == 3
        alpha = 0.55
        aa = 2.0^(2.0*alpha-1)
        aa = aa ^ (1/alpha)
        fi = aa * fi
        fi = fi .^ alpha
        for i = (no1+1):nbf5
            a = max(aa * n[i] * (1 - n[i]), 10^-15)
            for k = 1:nv
                dfi_dgamma[i, k] =
                    alpha * a^(alpha - 1) * aa * (1 - 2 * n[i]) * dn_dgamma[i, k]
            end
        end
    elseif ista == 4
        fi = sqrt.(fi)
        for i = (no1+1):nbf5
            a = max(fi[i], 10^-15)
            for k = 1:nv
                dfi_dgamma[i, k] =
                    1 / (2 * a) * (1 - 2 * n[i]) * dn_dgamma[i, k] * exp(-(n[i] - 0.5)^2)
                dfi_dgamma[i, k] +=
                    - fi[i] * 2 * (n[i] - 0.5) * dn_dgamma[i, k] * exp(-(n[i] - 0.5)^2)
            end
        end
        fi = fi .* exp.(- (n .- 0.5) .^ 2)
    end

    # Interpair Electron correlation #

    @tullio Dcj12r[i, j, k] := 2 * dn_dgamma[i, k] * n[j]
    @tullio Dck12r[i, j, k] := dn_dgamma[i, k] * n[j]

    if ista==0
        Dck12r_beta_cwo = view(Dck12r, (no1+1):nbeta, (nalpha+1):nbf5, 1:nv)
        Dck12r_cwo_beta = view(Dck12r, (nalpha+1):nbf5, (no1+1):nbeta, 1:nv)
        Dck12r_cwo_cwo = view(Dck12r, (nalpha+1):nbf5, (nalpha+1):nbf5, 1:nv)
        fi_beta = view(fi, (no1+1):nbeta)
        fi_cwo = view(fi, (nalpha+1):nbf5)
        dfi_beta = view(dfi_dgamma, (no1+1):nbeta, 1:nv)
        dfi_cwo = view(dfi_dgamma, (nalpha+1):nbf5, 1:nv)
        @tullio Dck12r_beta_cwo[i, j, k] += dfi_beta[i, k] * fi_cwo[j]
        @tullio Dck12r_cwo_beta[i, j, k] += fi_beta[j] * dfi_cwo[i, k]
        @tullio Dck12r_cwo_cwo[i, j, k] += fi_cwo[j] * dfi_cwo[i, k]

        if nsoc > 0
            Dck12r_beta_alpha = view(Dck12r, (no1+1):nbeta, (nbeta+1):nalpha, 1:nv)
            Dck12r_cwo_alpha = view(Dck12r, (nalpha+1):nbf5, (nbeta+1):nalpha, 1:nv)
            dfi_cwo = view(dfi_dgamma, (nalpha+1):nbf5, 1:nv)
            @tullio Dck12r_beta_alpha[i, j, k] += 0.5 * dfi_beta[i, k] * 0.5
            @tullio Dck12r_cwo_alpha[i, j, k] += dfi_cwo[i, k] * 0.5
        end

        if nsoc > 1
            Dck12r[(nbeta+1):nalpha, (nbeta+1):nalpha, 1:nv] .= 0.0
        end
    elseif ista>0
        @tullio Dck12r[i, j, k] := dn_dgamma[i, k] * n[j] + dfi_dgamma[i, k] * fi[j]
        Dck12r_beta_cwo = view(Dck12r, (no1+1):nbeta, (nalpha+1):nbf5, 1:nv)
        Dck12r_cwo_beta = view(Dck12r, (nalpha+1):nbf5, (no1+1):nbeta, 1:nv)
        Dck12r_cwo_cwo = view(Dck12r, (nalpha+1):nbf5, (nalpha+1):nbf5, 1:nv)
        if nsoc > 1
            Dck12r[(nbeta+1):nalpha, (nbeta+1):nalpha, 1:nv] .= 0.0
        end
    end

    n_d12_beta = view(n_d12, (no1+1):nbeta)
    n_d12_cwo = view(n_d12, (nalpha+1):nbf5)
    n_d_beta = view(n_d, (no1+1):nbeta)
    n_d_cwo = view(n_d, (nalpha+1):nbf5)
    dn_d12_beta = view(dn_d12_dgamma, (no1+1):nbeta, 1:nv)
    dn_d12_cwo = view(dn_d12_dgamma, (nalpha+1):nbf5, 1:nv)
    dn_d_beta = view(dn_d_dgamma, (no1+1):nbeta, 1:nv)
    dn_d_cwo = view(dn_d_dgamma, (nalpha+1):nbf5, 1:nv)
    @tullio Dck12r_beta_cwo[i, j, k] +=
        dn_d12_beta[i, k] * n_d12_cwo[j] - dn_d_beta[i, k] * n_d_cwo[j]
    @tullio Dck12r_cwo_beta[i, j, k] +=
        n_d12_beta[j] * dn_d12_cwo[i, k] - n_d_beta[j] * dn_d_cwo[i, k]
    @tullio Dck12r_cwo_cwo[i, j, k] +=
        -n_d12_cwo[j] * dn_d12_cwo[i, k] - n_d_cwo[j] * dn_d_cwo[i, k]

    # Intrapair Electron Correlation

    for l = 1:ndoc
        ldx = no1 + l

        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - l) + 1
        ul = ll + ndoc * (ncwo - 1)

        Dcj12r[ldx, ll:ndoc:ul, 1:nv] .= 0
        Dcj12r[ll:ndoc:ul, ldx, 1:nv] .= 0

        Dcj12r[ll:ndoc:ul, ll:ndoc:ul, 1:nv] .= 0

        a = max(n[ldx], 10^-15)
        b = n[ll:ndoc:ul]
        b[b .< 10^-15] .= 10^-15

        Dck12r_occ_cwo = view(Dck12r, ldx, ll:ndoc:ul, 1:nv)
        Dck12r_cwo_occ = view(Dck12r, ll:ndoc:ul, ldx, 1:nv)
        Dck12r_cwo_cwo = view(Dck12r, ll:ndoc:ul, ll:ndoc:ul, 1:nv)
        n_cwo = view(n, ll:ndoc:ul)
        n_occ = n[ldx]
        dn_dgamma_occ = view(dn_dgamma, ldx, 1:nv)
        dn_dgamma_cwo = view(dn_dgamma, ll:ndoc:ul, 1:nv)
        @tullio Dck12r_occ_cwo[i, j] =
            1 / 2 * 1 / sqrt(a) * dn_dgamma_occ[j] * sqrt(n_cwo[i])
        @tullio Dck12r_cwo_occ[i, j] =
            1 / 2 * 1 / sqrt(b[i]) * dn_dgamma_cwo[i, j] * sqrt(n_occ)
        @tullio Dck12r_cwo_cwo[i, j, k] =
            -1 / 2 * 1 / sqrt(b[i]) * dn_dgamma_cwo[i, k] * sqrt(n_cwo[j])

    end

    return Dcj12r, Dck12r

end

function CJCKD9(n, no1, ndoc, nsoc, nbeta, nalpha, ndns, ncwo, MSpin, h_cut)

    nbf5 = size(n)[1]

    #h_cut = p.h_cut#0.02*sqrt(2.0)
    n_d = zeros(nbf5)

    for i = 1:ndoc
        idx = no1 + i
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - i) + 1
        ul = ll + ndoc * (ncwo - 1)

        h = 1.0 - n[idx]
        coc = h / h_cut
        arg = -coc^2
        F = exp(arg)  # ! Hd/Hole
        n_d[idx] = n[idx] * F
        n_d[ll:ndoc:ul] = n[ll:ndoc:ul] * F  # ROd = RO*Hd/Hole

    end

    n_d12 = sqrt.(n_d)
    fi = n .* (1 .- n)
    fi[fi .<= 0] .= 0
    fi = sqrt.(fi)

    # Interpair Electron correlation #

    @tullio cj12[i, j] := 2 * n[i] * n[j]

    #@tullio ck12[i,j] := n[i]*n[j]
    @tullio ck12[i, j] := n[i] * n[j] + fi[i] * fi[j]
    ck12_beta_cwo = view(ck12, (no1+1):nbeta, (nalpha+1):nbf5)
    ck12_cwo_beta = view(ck12, (nalpha+1):nbf5, (no1+1):nbeta)
    ck12_cwo_cwo = view(ck12, (nalpha+1):nbf5, (nalpha+1):nbf5)
    #fi_beta = view(fi,no1+1:nbeta)
    #fi_cwo = view(fi,nalpha+1:nbf5)
    #@tullio ck12_beta_cwo[i,j] += fi_beta[i]*fi_cwo[j]
    #@tullio ck12_cwo_beta[i,j] += fi_cwo[i]*fi_beta[j]
    #@tullio ck12_cwo_cwo[i,j] += fi_cwo[i]*fi_cwo[j]

    # Intrapair Electron Correlation

    #if(MSpin==0 && nsoc>0)
    #    ck12[no1+1:nbeta,nbeta+1:nalpha] .+= 0.5*fi[no1+1:nbeta]*0.5
    #    ck12[nbeta+1:nalpha,no1+1:nbeta] .+= 0.5*0.5*fi[no1+1:nbeta]'
    #    ck12[nbeta+1:nalpha,nalpha+1:nbf5] .+= 0.5*fi[nalpha+1:nbf5]'
    #    ck12[nalpha+1:nbf5,nbeta+1:nalpha] .+= fi[nalpha+1:nbf5]*0.5
    #end

    if nsoc > 1
        ck12[(nbeta+1):nalpha, (nbeta+1):nalpha] .= 0.5
    end

    n_d12_beta = view(n_d12, (no1+1):nbeta)
    n_d12_cwo = view(n_d12, (nalpha+1):nbf5)
    n_d_beta = view(n_d, (no1+1):nbeta)
    n_d_cwo = view(n_d, (nalpha+1):nbf5)
    @tullio ck12_beta_cwo[i, j] += n_d12_beta[i] * n_d12_cwo[j] - n_d_beta[i] * n_d_cwo[j]
    @tullio ck12_cwo_beta[i, j] += n_d12_cwo[i] * n_d12_beta[j] - n_d_cwo[i] * n_d_beta[j]
    @tullio ck12_cwo_cwo[i, j] += -n_d12_cwo[i] * n_d12_cwo[j] - n_d_cwo[i] * n_d_cwo[j]

    for l = 1:ndoc
        ldx = no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - l) + 1
        ul = ll + ndoc * (ncwo - 1)

        cj12[ldx, ll:ndoc:ul] .= 0
        cj12[ll:ndoc:ul, ldx] .= 0

        cj12[ll:ndoc:ul, ll:ndoc:ul] .= 0

        ck12[ldx, ll:ndoc:ul] .= sqrt.(n[ldx] * n[ll:ndoc:ul])
        ck12[ll:ndoc:ul, ldx] .= sqrt.(n[ldx] * n[ll:ndoc:ul])

        ck12_ww = view(ck12, ll:ndoc:ul, ll:ndoc:ul)
        n_ww = view(n, ll:ndoc:ul)
        @tullio ck12_ww[i, j] = -sqrt(n_ww[i] * n_ww[j])
    end

    return cj12, ck12

end

function der_CJCKD9(
    n,
    dn_dgamma,
    no1,
    ndoc,
    nalpha,
    nbeta,
    nv,
    nbf5,
    ndns,
    ncwo,
    MSpin,
    nsoc,
    h_cut,
)

    nbf5 = size(n)[1]

    #h_cut = p.h_cut#0.02*sqrt(2.0)
    n_d = zeros(nbf5)
    dn_d_dgamma = zeros(nbf5, nv)
    dn_d12_dgamma = zeros(nbf5, nv)

    for i = 1:ndoc
        idx = no1 + i
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - i) + 1
        ul = ll + ndoc * (ncwo - 1)

        h_idx = 1.0 - n[idx]
        coc = h_idx / h_cut
        arg = -coc^2
        F_idx = exp(arg)                # Hd/Hole
        n_d[idx] = n[idx] * F_idx
        n_d[ll:ndoc:ul] = n[ll:ndoc:ul] * F_idx      # n_d = RO*Hd/Hole
        dn_d_dgamma[idx, 1:nv] .=
            F_idx * dn_dgamma[idx, 1:nv] .* (1 - n[idx] * (-2 * coc / h_cut))
        dn_d_dgamma[ll:ndoc:ul, 1:nv] .= F_idx * dn_dgamma[ll:ndoc:ul, 1:nv]
        dn_d_dgamma_cwo_v = view(dn_d_dgamma, ll:ndoc:ul, 1:nv)
        dn_dgamma_v = view(dn_dgamma, idx, 1:nv)
        n_cwo = view(n, ll:ndoc:ul)
        @tullio dn_d_dgamma_cwo_v[i, j] +=
            F_idx * (2 * coc / h_cut) * n_cwo[i] * dn_dgamma_v[j]
    end

    n_d12 = sqrt.(n_d)
    dn_d12_dgamma = 0.5 * dn_d_dgamma ./ max.(n_d12, 10^-15)

    fi = n .* (1 .- n)
    fi[fi .<= 0] .= 0
    fi = sqrt.(fi)

    dfi_dgamma = zeros(nbf5, nv)
    for i = (no1+1):nbf5
        a = max(fi[i], 10^-15)
        for k = 1:nv
            dfi_dgamma[i, k] = 1 / (2 * a) * (1 - 2 * n[i]) * dn_dgamma[i, k]
        end
    end


    # Interpair Electron correlation #

    @tullio Dcj12r[i, j, k] := 2 * dn_dgamma[i, k] * n[j]

    #@tullio Dck12r[i,j,k] := dn_dgamma[i,k]*n[j]
    @tullio Dck12r[i, j, k] := dn_dgamma[i, k] * n[j] + dfi_dgamma[i, k] * fi[j]
    Dck12r_beta_cwo = view(Dck12r, (no1+1):nbeta, (nalpha+1):nbf5, 1:nv)
    Dck12r_cwo_beta = view(Dck12r, (nalpha+1):nbf5, (no1+1):nbeta, 1:nv)
    Dck12r_cwo_cwo = view(Dck12r, (nalpha+1):nbf5, (nalpha+1):nbf5, 1:nv)
    #fi_beta = view(fi,no1+1:nbeta)
    #fi_cwo = view(fi,nalpha+1:nbf5)
    #dfi_beta = view(dfi_dgamma,no1+1:nbeta,1:nv)
    #dfi_cwo = view(dfi_dgamma,nalpha+1:nbf5,1:nv)
    #@tullio Dck12r_beta_cwo[i,j,k] += dfi_beta[i,k]*fi_cwo[j]
    #@tullio Dck12r_cwo_beta[i,j,k] += fi_beta[j]*dfi_cwo[i,k]
    #@tullio Dck12r_cwo_cwo[i,j,k] += fi_cwo[j]*dfi_cwo[i,k]

    #if(MSpin==0 && nsoc>0)
    #    Dck12r_beta_alpha = view(Dck12r,no1+1:nbeta,nbeta+1:nalpha,1:nv)
    #    Dck12r_cwo_alpha = view(Dck12r,nalpha+1:nbf5,nbeta+1:nalpha,1:nv)
    #    dfi_cwo = view(dfi_dgamma,nalpha+1:nbf5,1:nv)
    #	@tullio Dck12r_beta_alpha[i,j,k] += 0.5*dfi_beta[i,k]*0.5
    #	@tullio Dck12r_cwo_alpha[i,j,k] += dfi_cwo[i,k]*0.5
    #end

    if nsoc > 1
        Dck12r[(nbeta+1):nalpha, (nbeta+1):nalpha, 1:nv] .= 0.0
    end

    n_d12_beta = view(n_d12, (no1+1):nbeta)
    n_d12_cwo = view(n_d12, (nalpha+1):nbf5)
    n_d_beta = view(n_d, (no1+1):nbeta)
    n_d_cwo = view(n_d, (nalpha+1):nbf5)
    dn_d12_beta = view(dn_d12_dgamma, (no1+1):nbeta, 1:nv)
    dn_d12_cwo = view(dn_d12_dgamma, (nalpha+1):nbf5, 1:nv)
    dn_d_beta = view(dn_d_dgamma, (no1+1):nbeta, 1:nv)
    dn_d_cwo = view(dn_d_dgamma, (nalpha+1):nbf5, 1:nv)
    @tullio Dck12r_beta_cwo[i, j, k] +=
        dn_d12_beta[i, k] * n_d12_cwo[j] - dn_d_beta[i, k] * n_d_cwo[j]
    @tullio Dck12r_cwo_beta[i, j, k] +=
        n_d12_beta[j] * dn_d12_cwo[i, k] - n_d_beta[j] * dn_d_cwo[i, k]
    @tullio Dck12r_cwo_cwo[i, j, k] +=
        -n_d12_cwo[j] * dn_d12_cwo[i, k] - n_d_cwo[j] * dn_d_cwo[i, k]

    # Intrapair Electron Correlation

    for l = 1:ndoc
        ldx = no1 + l

        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + (ndoc - l) + 1
        ul = ll + ndoc * (ncwo - 1)

        Dcj12r[ldx, ll:ndoc:ul, 1:nv] .= 0
        Dcj12r[ll:ndoc:ul, ldx, 1:nv] .= 0

        Dcj12r[ll:ndoc:ul, ll:ndoc:ul, 1:nv] .= 0

        a = max(n[ldx], 10^-15)
        b = n[ll:ndoc:ul]
        b[b .< 10^-15] .= 10^-15

        Dck12r_occ_cwo = view(Dck12r, ldx, ll:ndoc:ul, 1:nv)
        Dck12r_cwo_occ = view(Dck12r, ll:ndoc:ul, ldx, 1:nv)
        Dck12r_cwo_cwo = view(Dck12r, ll:ndoc:ul, ll:ndoc:ul, 1:nv)
        n_cwo = view(n, ll:ndoc:ul)
        n_occ = n[ldx]
        dn_dgamma_occ = view(dn_dgamma, ldx, 1:nv)
        dn_dgamma_cwo = view(dn_dgamma, ll:ndoc:ul, 1:nv)
        @tullio Dck12r_occ_cwo[i, j] =
            1 / 2 * 1 / sqrt(a) * dn_dgamma_occ[j] * sqrt(n_cwo[i])
        @tullio Dck12r_cwo_occ[i, j] =
            1 / 2 * 1 / sqrt(b[i]) * dn_dgamma_cwo[i, j] * sqrt(n_occ)
        @tullio Dck12r_cwo_cwo[i, j, k] =
            -1 / 2 * 1 / sqrt(b[i]) * dn_dgamma_cwo[i, k] * sqrt(n_cwo[j])

    end

    return Dcj12r, Dck12r

end

function ocupacion(gamma, no1, ndoc, nalpha, nv, nbf5, ndns, ncwo, HighSpin, occ_method)
    if occ_method == "Trigonometric"
        n, dn_dgamma = ocupacion_trigonometric(
            gamma,
            no1,
            ndoc,
            nalpha,
            nv,
            nbf5,
            ndns,
            ncwo,
            HighSpin,
        )
    elseif occ_method == "Softmax"
        n, dn_dgamma =
            ocupacion_softmax(gamma, no1, ndoc, nalpha, nv, nbf5, ndns, ncwo, HighSpin)
    elseif occ_method == "EBI"
        n, dn_dgamma =
            ocupacion_ebi(gamma, no1, ndoc, nalpha, nv, nbf5, ndns, ncwo, HighSpin)
    end
    return n, dn_dgamma
end

function ocupacion_trigonometric(gamma, no1, ndoc, nalpha, nv, nbf5, ndns, ncwo, HighSpin)

    n = zeros(nbf5)
    dni_dgammai = zeros(nbf5)
    dn_dgamma = zeros(nbf5, nv)

    n[1:no1] .= 1

    n[(no1+1):(no1+ndoc)] .= 1 / 2 * (1 .+ (cos.(gamma[1:ndoc])) .^ 2)
    dni_dgammai[(no1+1):(no1+ndoc)] .= -1 / 2 * sin.(2 * gamma[1:ndoc])

    if !HighSpin
        n[(no1+ndoc+1):(no1+ndns)] .= 0.5
    elseif HighSpin
        n[(no1+ndoc+1):(no1+ndns)] .= 1.0
    end

    dn_dgamma = zeros(nbf5, nv)
    h = 1 .- n

    for i = 1:ndoc
        ll_n = no1 + ndns + (ndoc - i) + 1
        ul_n = ll_n + ndoc * (ncwo - 1)
        n_pi = @view n[ll_n:ndoc:ul_n]

        ll_gamma = ndoc + (ndoc - i) + 1
        ul_gamma = ll_gamma + ndoc * (ncwo - 2)
        gamma_pi = @view gamma[ll_gamma:ndoc:ul_gamma]

        n_pi .= h[no1+i]
        for kw = 1:(ncwo-1)
            n_pi[kw] *= sin(gamma_pi[kw])^2
            n_pi[(kw+1):end] .*= cos(gamma_pi[kw])^2
        end

        # dn_g/dgamma_g
        dn_dgamma[no1+i, i] = dni_dgammai[no1+i]

        # dn_pi/dgamma_g
        dn_pi_dgamma_g = @view dn_dgamma[ll_n:ndoc:ul_n, i]
        dn_pi_dgamma_g .= -dni_dgammai[no1+i]
        for kw = 1:(ncwo-1)
            dn_pi_dgamma_g[kw] *= sin(gamma_pi[kw])^2
            dn_pi_dgamma_g[(kw+1):end] .*= cos(gamma_pi[kw])^2
        end

        # dn_pi/dgamma_pj (j<i)
        dn_pi_dgamma_pj = @view dn_dgamma[ll_n:ndoc:ul_n, ll_gamma:ndoc:ul_gamma]
        for jw = 1:(ncwo-1)
            dn_pi_dgamma_pj[(jw+1):end, jw] .= n[no1+i] - 1
            for kw = 1:(jw-1)
                dn_pi_dgamma_pj[(jw+1):end, jw] .*= cos(gamma_pi[kw])^2
            end
            dn_pi_dgamma_pj[(jw+1):end, jw] .*= sin(2 * gamma_pi[jw])
            for kw = (jw+1):(ncwo-1)
                dn_pi_dgamma_pj[kw, jw] *= sin(gamma_pi[kw])^2
                dn_pi_dgamma_pj[(kw+1):end, jw] .*= cos(gamma_pi[kw])^2
            end
        end

        # dn_pi/dgamma_pi
        for jw = 1:(ncwo-1)
            dn_pi_dgamma_pj[jw, jw] = 1 - n[no1+i]
        end
        for kw = 1:(ncwo-1)
            dn_pi_dgamma_pj[kw, kw] *= sin(2 * gamma_pi[kw])
            for lw = (kw+1):(ncwo-1)
                dn_pi_dgamma_pj[lw, lw] *= cos(gamma_pi[kw])^2
            end
        end
    end

    return n, dn_dgamma
end

function ocupacion_softmax(gamma, no1, ndoc, nalpha, nv, nbf5, ndns, ncwo, HighSpin)

    n = Vector{Float64}(undef, nbf5)
    dn_dgamma = zeros(nbf5, nv)

    n[1:no1] .= 1                                              # [1,no1]
    if !HighSpin
        n[(no1+ndoc+1):(no1+ndns)] .= 0.5   # (no1+ndoc,no1+ndns]
    elseif HighSpin
        n[(no1+ndoc+1):(no1+ndns)] .= 1.0   # (no1+ndoc,no1+ndns]
    end

    exp_gamma = exp.(gamma)

    for i = 1:ndoc

        ll = no1 + ndns + (ndoc - i) + 1
        ul = ll + ndoc * (ncwo - 1)
        n_pi = @view n[ll:ndoc:ul]

        llg = (ndoc - i) + 1
        ulg = llg + ndoc * (ncwo - 1)
        exp_gamma_pi = @view exp_gamma[llg:ndoc:ulg]

        den = 1 + sum(exp_gamma_pi)

        n[no1+i] = 1 / den
        n_pi .= exp_gamma_pi / den

        dn_g_dgamma_pi = @view dn_dgamma[no1 + i, llg:ndoc:ulg]
        dn_g_dgamma_pi .= - n[i] * n_pi

        dn_pi_dgamma_pi = @view dn_dgamma[ll:ndoc:ul, llg:ndoc:ulg]
        dn_pi_dgamma_pi .= - n_pi * n_pi'

        dn_pi_dgamma_pi = @view dn_dgamma[ll:ndoc:ul, llg:ndoc:ulg]
        dn_pi_dgamma_pi[diagind(dn_pi_dgamma_pi)] += n_pi

    end

    return n, dn_dgamma
end

function ocupacion_ebi(x, no1, ndoc, nalpha, nv, nbf5, ndns, ncwo, HighSpin)
    """Transform gammas to n according to the ebi
    parameterization of the occupation numbers

    n_p = (1 + erf(x_p + mu)) / 2

    dn_p/dx_p = erf'(x_p + mu) * (1 + dmu/dx_p) / 2
    dn_p/dx_q = erf'(x_p + mu) * dmu/dx_q / 2

    dmu/dx_p = - erf'(x_p + mu)/(sum_q erf'(x_q + mu)) 
    """

    nv = size(x)[1]

    mu = [0.0]
    N = ndoc
    res = optimize(
        mu -> deviation_ebi(mu, x, N),
        mu -> D_deviation_ebi(mu, x, N),
        mu,
        ConjugateGradient(),
        inplace = false,
    )
    #res = optimize(mu->deviation_ebi(mu,x,N), mu->D_deviation_ebi(mu,x,N), mu, SimulatedAnnealing(), inplace=false)
    mu = res.minimizer

    # n_p
    n = zeros(nv)
    for (i, xi) in enumerate(x)
        n[i] = (1.0 + erf(xi + mu[1])) / 2
    end

    #println("sum n")
    #println((n))
    #println(sum(n))

    # dmu/dx_p
    den = 0
    for i = 1:nv
        den += exp(-(x[i] + mu[1])^2)
    end
    D_mu_x = zeros(nv)
    for i = 1:nv
        D_mu_x[i] = -exp(-(x[i] + mu[1])^2) / den
    end

    # dn_p/dx_
    D_n_x = zeros(nv, nv)
    for i = 1:nv
        for j = 1:nv
            if i == j
                D_n_x[i, j] = 1 / sqrt(pi) * exp(-(x[i] + mu[1])^2) * (1 + D_mu_x[j])
            else
                D_n_x[i, j] = 1 / sqrt(pi) * exp(-(x[i] + mu[1])^2) * D_mu_x[j]
            end
        end
    end

    return n, D_n_x
end

function deviation_ebi(mu, x, N)
    """Loss function used in ebi: Squared deviation from N/2 of the
    occupation numbers

    f = (N - sum_i n_i)^2

    """

    nv = size(x)[1]

    n = zeros(nv)
    for i = 1:nv
        n[i] += (1 + erf(x[i] + mu[1])) / 2
    end
    sum_n = sum(n)

    dev = (N - sum_n)^2

    return dev
end

function D_deviation_ebi(mu, x, N)
    """Derivative from gamma of the loss function used in ebi

    df/dmu = 2 * (N- sum_p n_p) * (-sum_p dn_p/dmu)

    """

    nv = size(x)[1]

    # n_p
    sum_n = 0
    for i = 1:nv
        sum_n += (1 + erf(x[i] + mu[1])) / 2
    end

    # dn_p/dmu = 1/sqrt(pi) * exp(-(x_p + mu)^2)
    sum_dn = 0.0
    for i = 1:nv
        sum_dn += 1 / sqrt(pi) * exp(-(x[i] + mu[1])^2)
    end

    D_dev = -2 * (N - sum_n) * sum_dn# + 0.1*(N-sum_n)

    return D_dev
end

function calce(n, cj12, ck12, J_MO, K_MO, H_core, p)

    E = 0

    n_beta = view(n, 1:p.nbeta)
    n_alpha = view(n, (p.nbeta+1):p.nalpha)
    n_nbf5 = view(n, (p.nalpha+1):p.nbf5)
    H_beta = view(H_core, 1:p.nbeta)
    H_alpha = view(H_core, (p.nbeta+1):p.nalpha)
    H_nbf5 = view(H_core, (p.nalpha+1):p.nbf5)
    J_MO_beta = view(J_MO, 1:p.nbeta, 1:p.nbeta)
    J_MO_nbf5 = view(J_MO, (p.nalpha+1):p.nbf5, (p.nalpha+1):p.nbf5)

    # 2H + J
    @tullio E += n_beta[i] * (2 * H_beta[i] + J_MO_beta[i, i])
    if (p.nsoc > 0)
        @tullio E += n_alpha[i] * 2 * H_alpha[i]
    end
    @tullio E += n_nbf5[i] * (2 * H_nbf5[i] + J_MO_nbf5[i, i])

    #C^J JMO
    @tullio E += cj12[i, j] * J_MO[j, i]

    #C^K KMO
    @tullio E += -ck12[i, j] * K_MO[j, i]

    return E
end

function calcocce(gamma, J_MO, K_MO, H_core, p)

    n, dn_dgamma = ocupacion(
        gamma,
        p.no1,
        p.ndoc,
        p.nalpha,
        p.nv,
        p.nbf5,
        p.ndns,
        p.ncwo,
        p.HighSpin,
        p.occ_method,
    )
    cj12, ck12 = PNOFi_selector(n, p)

    E = calce(n, cj12, ck12, J_MO, K_MO, H_core, p)

    return E

end

function calcoccg(gamma, J_MO, K_MO, H_core, p)

    grad = zeros(p.nv)
    n, dn_dgamma = ocupacion(
        gamma,
        p.no1,
        p.ndoc,
        p.nalpha,
        p.nv,
        p.nbf5,
        p.ndns,
        p.ncwo,
        p.HighSpin,
        p.occ_method,
    )
    Dcj12r, Dck12r = der_PNOFi_selector(n, dn_dgamma, p)

    # dn_dgamma (2H+J)
    dn_dgamma_beta = view(dn_dgamma, (p.no1+1):p.nbeta, 1:p.nv)
    dn_dgamma_nbf5 = view(dn_dgamma, (p.nalpha+1):p.nbf5, 1:p.nv)
    H_core_beta = view(H_core, (p.no1+1):p.nbeta)
    H_core_nbf5 = view(H_core, (p.nalpha+1):p.nbf5)
    J_MO_beta = view(J_MO, (p.no1+1):p.nbeta, (p.no1+1):p.nbeta)
    J_MO_nbf5 = view(J_MO, (p.nalpha+1):p.nbf5, (p.nalpha+1):p.nbf5)
    @tullio grad[k] += dn_dgamma_beta[i, k] * 2 * H_core_beta[i]# + J_MO_beta[i,i])
    @tullio grad[k] += dn_dgamma_beta[i, k] * J_MO_beta[i, i]
    #grad += np.einsum('ik,i->k',dn_dgamma[p.no1:p.nbeta,:p.nv],2*H_core[p.no1:p.nbeta]+np.diagonal(J_MO)[p.no1:p.nbeta],optimize=True) # [0,Nbeta]
    @tullio grad[k] += dn_dgamma_nbf5[i, k] * 2 * H_core_nbf5[i]# + J_MO_nbf5[i,i])
    @tullio grad[k] += dn_dgamma_nbf5[i, k] * J_MO_nbf5[i, i]

    # 2 dCJ_dgamma J_MO
    Dcj12r_beta = view(Dcj12r, (p.no1+1):p.nbeta, 1:p.nbf5, 1:p.nv)
    Dcj12r_nbf5 = view(Dcj12r, (p.nalpha+1):p.nbf5, 1:p.nbf5, 1:p.nv)
    J_MO_beta = view(J_MO, 1:p.nbf5, (p.no1+1):p.nbeta)
    J_MO_nbf5 = view(J_MO, 1:p.nbf5, (p.nalpha+1):p.nbf5)
    @tullio grad[k] += 2 * Dcj12r_beta[i, j, k] * J_MO_beta[j, i]
    #grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.no1:p.nbeta,:p.nbf5,:p.nv],J_MO[:p.nbf5,p.no1:p.nbeta],optimize=True)

    @tullio grad[k] += 2 * Dcj12r_nbf5[i, j, k] * J_MO_nbf5[j, i]
    #grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.nalpha:p.nbf5,:p.nbf5,:p.nv],J_MO[:p.nbf5,p.nalpha:p.nbf5],optimize=True)

    # -2 dCK_dgamma K_MO
    Dck12r_beta = view(Dck12r, (p.no1+1):p.nbeta, 1:p.nbf5, 1:p.nv)
    Dck12r_nbf5 = view(Dck12r, (p.nalpha+1):p.nbf5, 1:p.nbf5, 1:p.nv)
    K_MO_beta = view(K_MO, 1:p.nbf5, (p.no1+1):p.nbeta)
    K_MO_nbf5 = view(K_MO, 1:p.nbf5, (p.nalpha+1):p.nbf5)
    @tullio grad[k] += -2 * Dck12r_beta[i, j, k] * K_MO_beta[j, i]
    #grad -= 2*np.einsum('ijk,ji->k',Dck12r[p.no1:p.nbeta,:p.nbf5,:p.nv],K_MO[:p.nbf5,p.no1:p.nbeta],optimize=True)

    @tullio grad[k] += -2 * Dck12r_nbf5[i, j, k] * K_MO_nbf5[j, i]
    #grad -= 2*np.einsum('ijk,ji->k',Dck12r[p.nalpha:p.nbf5,:p.nbf5,:p.nv],K_MO[:p.nbf5,p.nalpha:p.nbf5],optimize=True)

    return grad
end

function calcorbe(y, n, cj12, ck12, C, H, I, b_mnl, p)

    Cnew = rotate_orbital(y, C, p)
    J_MO, K_MO, H_core = computeJKH_MO(Cnew, H, I, b_mnl, p)
    E = calce(n, cj12, ck12, J_MO, K_MO, H_core, p)

    return E

end

function calcorbg(y, n, cj12, ck12, C, H, I, pa)

    Cnew = rotate_orbital(y, C, pa)

    elag, Hmat = compute_Lagrange2(Cnew, n, H, I, cj12, ck12, pa.nalpha, pa.nbeta)

    grad = 4 * elag - 4 * elag'
    grads = zeros(pa.nvar)
    n = 1
    for i = 1:pa.nbf5
        for j = (i+1):pa.nbf
            grads[n] = grad[i, j]
            n += 1
        end
    end

    return grads

end

function calcorbeg(F, G, y, n, cj12, ck12, C, H, I, pa)

    Cnew = rotate_orbital(y, C, pa)

    elag, Hmat = compute_Lagrange2(Cnew, n, H, I, cj12, ck12, pa.nalpha, pa.nbeta)

    if G !== nothing
        grad = 4 * elag - 4 * elag'
        grads = zeros(pa.nvar)
        nn = 1
        for i = 1:pa.nbf5
            for j = (i+1):pa.nbf
                grads[nn] = grad[i, j]
                nn += 1
            end
        end
        G .= grads
    end
    if F !== nothing
        E = computeE_elec(Hmat, n, elag, pa)
        return E
    end

end

function calccombe(x, C, H, I, b_mnl, p)

    y = x[1:p.nvar]
    gamma = x[(p.nvar+1):end]

    Cnew = rotate_orbital(y, C, p)

    n, dn_dgamma = ocupacion(
        gamma,
        p.no1,
        p.ndoc,
        p.nalpha,
        p.nv,
        p.nbf5,
        p.ndns,
        p.ncwo,
        p.HighSpin,
        p.occ_method,
    )
    cj12, ck12 = PNOFi_selector(n, p)

    J_MO, K_MO, H_core = computeJKH_MO(Cnew, H, I, b_mnl, p)

    E = calce(n, cj12, ck12, J_MO, K_MO, H_core, p)

    return E
end

function calccombg(x, C, H, I, b_mnl, p)

    y = x[1:p.nvar]
    gamma = x[(p.nvar+1):end]

    Cnew = rotate_orbital(y, C, p)

    n, dn_dgamma = ocupacion(
        gamma,
        p.no1,
        p.ndoc,
        p.nalpha,
        p.nv,
        p.nbf5,
        p.ndns,
        p.ncwo,
        p.HighSpin,
        p.occ_method,
    )
    cj12, ck12 = PNOFi_selector(n, p)

    J_MO, K_MO, H_core = computeJKH_MO(Cnew, H, I, b_mnl, p)

    grad = zeros(p.nvar + p.nv)

    grad[1:p.nvar] = calcorbg(y, n, cj12, ck12, C, H, I, b_mnl, p)
    grad[(p.nvar+1):end] = calcoccg(gamma, J_MO, K_MO, H_core, p)

    return grad
end

function calccombeg(F, G, x, C, H, I, b_mnl, p)

    y = x[1:p.nvar]
    gamma = x[(p.nvar+1):end]

    Cnew = rotate_orbital(y, C, p)

    n, dn_dgamma = ocupacion(
        gamma,
        p.no1,
        p.ndoc,
        p.nalpha,
        p.nv,
        p.nbf5,
        p.ndns,
        p.ncwo,
        p.HighSpin,
        p.occ_method,
    )
    cj12, ck12 = PNOFi_selector(n, p)

    J_MO, K_MO, H_core = computeJKH_MO(Cnew, H, I, b_mnl, p)

    if G !== nothing
        grad = zeros(p.nvar + p.nv)

        grad[(p.nvar+1):end] = calcoccg(gamma, J_MO, K_MO, H_core, p)
        if F !== nothing
            gv = view(grad, 1:p.nvar)
            E = calcorbeg(F, gv, y, n, cj12, ck12, C, H, I, b_mnl, p)
            G .= grad
            return E
        else
            grad[1:p.nvar] = calcorbg(y, n, cj12, ck12, C, H, I, b_mnl, p)
            G .= grad
        end
    end
    if F !== nothing
        E = calcocce(gamma, J_MO, K_MO, H_core, p)
        return E
    end

end

function compute_2RDM(pp, n)

    # PNOF5
    # Interpair Electron correlation
    Id = 1.0 * Matrix(I, pp.nbf5, pp.nbf5)

    inter = n .* n'
    intra = zeros(pp.nbf5, pp.nbf5)

    # Intrapair Electron Correlation
    for l = 1:pp.ndoc
        ldx = pp.no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = pp.no1 + pp.ndns + (pp.ndoc - l) + 1
        ul = ll + pp.ndoc * (pp.ncwo - 1)

        inter[ldx, ldx] = 0
        inter[ldx, ll:pp.ndoc:ul] .= 0
        inter[ll:pp.ndoc:ul, ldx] .= 0
        inter[ll:pp.ndoc:ul, ll:pp.ndoc:ul] .= 0

        intra[ldx, ldx] = sqrt(n[ldx] * n[ldx])
        intra[ldx, ll:pp.ndoc:ul] .= -sqrt.(n[ldx] * n[ll:pp.ndoc:ul])
        intra[ll:pp.ndoc:ul, ldx] .= -sqrt.(n[ldx] * n[ll:pp.ndoc:ul])

        intra_ww = view(intra, ll:pp.ndoc:ul, ll:pp.ndoc:ul)
        n_ww = view(n, ll:pp.ndoc:ul)
        @tullio intra_ww[i, j] = sqrt(n_ww[i] * n_ww[j])
    end

    for i = (pp.nbeta+1):pp.nalpha
        inter[i, i] = 0
    end

    @tullio Daa[p, q, r, t] := inter[p, q] * Id[p, r] * Id[q, t]
    @tullio Daa[p, q, r, t] += -inter[p, q] * Id[p, t] * Id[q, r]
    @tullio Dab[p, q, r, t] := inter[p, q] * Id[p, r] * Id[q, t]
    @tullio Dab[p, q, r, t] += intra[p, r] * Id[p, q] * Id[r, t]

    # PNOF7
    if (pp.ipnof == 7 || pp.ipnof == 8)
        fi = n .* (1 .- n)
        fi[fi .<= 0] .= 0
        if (pp.ista == 2)
            fi = 0.9 * sqrt.(fi)
        else
            fi = sqrt.(fi)
        end
        Pi_s = fi .* fi'
        # Intrapair Electron Correlation
        for l = 1:pp.ndoc
            ldx = pp.no1 + l
            # inicio y fin de los orbitales acoplados a los fuertemente ocupados
            ll = pp.no1 + pp.ndns + (pp.ndoc - l) + 1
            ul = ll + pp.ndoc * (pp.ncwo - 1)

            Pi_s[ldx, ldx] = 0
            Pi_s[ldx, ll:pp.ndoc:ul] .= 0
            Pi_s[ll:pp.ndoc:ul, ldx] .= 0
            Pi_s[ll:pp.ndoc:ul, ll:pp.ndoc:ul] .= 0
        end
        for i = (pp.nbeta+1):pp.nalpha
            Pi_s[i, i] = 0
        end

        inter2 = 0 * Pi_s
        inter2[(pp.nbeta+1):pp.nalpha, (pp.nbeta+1):pp.nalpha] =
            Pi_s[(pp.nbeta+1):pp.nalpha, (pp.nbeta+1):pp.nalpha]
        Pi_s[(pp.nbeta+1):pp.nalpha, (pp.nbeta+1):pp.nalpha] .= 0

        @tullio Dab[p, q, r, t] += -inter2[p, q] * Id[p, t] * Id[q, r]

        if (pp.ipnof == 8 && pp.ista==0)
            Pi_s[1:pp.nbeta, 1:pp.nbeta] .= 0
            Pi_s[1:pp.nbeta, (pp.nbeta+1):pp.nalpha] *= 0.5
            Pi_s[(pp.nbeta+1):pp.nalpha, 1:pp.nbeta] *= 0.5
        end

        @tullio Dab[p, q, r, t] += -Pi_s[p, r] * Id[p, q] * Id[r, t]

        if (pp.ipnof == 8)

            h_cut = pp.h_cut
            n_d = zeros(size(n)[1])

            for i = 1:pp.ndoc
                idx = pp.no1 + i
                # inicio y fin de los orbitales acoplados a los fuertemente ocupados
                ll = pp.no1 + pp.ndns + (pp.ndoc - i) + 1
                ul = ll + pp.ndoc * (pp.ncwo - 1)

                h = 1.0 - n[idx]
                coc = h / h_cut
                arg = -coc^2
                F = exp(arg)  # ! Hd/Hole
                n_d[idx] = n[idx] * F
                n_d[ll:pp.ndoc:ul] = n[ll:pp.ndoc:ul] * F  # ROd = RO*Hd/Hole
            end

            n_d12 = sqrt.(n_d)

            inter = n_d12 .* n_d12' - n_d .* n_d'
            inter2 = n_d12 .* n_d12' + n_d .* n_d'
            # Intrapair Electron Correlation
            for l = 1:pp.ndoc
                ldx = pp.no1 + l
                # inicio y fin de los orbitales acoplados a los fuertemente ocupados
                ll = pp.no1 + pp.ndns + (pp.ndoc - l) + 1
                ul = ll + pp.ndoc * (pp.ncwo - 1)

                inter[ldx, ldx] = 0
                inter[ldx, ll:pp.ndoc:ul] .= 0
                inter[ll:pp.ndoc:ul, ldx] .= 0
                inter[ll:pp.ndoc:ul, ll:pp.ndoc:ul] .= 0
                inter2[ldx, ldx] = 0
                inter2[ldx, ll:pp.ndoc:ul] .= 0
                inter2[ll:pp.ndoc:ul, ldx] .= 0
                inter2[ll:pp.ndoc:ul, ll:pp.ndoc:ul] .= 0
            end

            inter[(pp.nbeta+1):end, (pp.nbeta+1):end] .= 0
            inter[1:pp.nalpha, 1:pp.nalpha] .= 0
            inter2[1:pp.nalpha, :] .= 0
            inter2[:, 1:pp.nalpha] .= 0

            Pi_d = inter - inter2

            @tullio Dab[p, q, r, t] += -Pi_d[p, r] * Id[p, q] * Id[r, t]
        end
    end

    return Daa, Dab

end
