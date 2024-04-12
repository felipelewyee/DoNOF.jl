function PNOFi_selector(n,p)
    if p.ipnof==5
        cj12,ck12 = CJCKD5(n,p.no1,p.ndoc,p.nsoc,p.nbeta,p.nalpha,p.ndns,p.ncwo,p.MSpin)
    elseif p.ipnof==7
        cj12,ck12 = CJCKD7(n,p.ista,p.no1,p.ndoc,p.nsoc,p.nbeta,p.nalpha,p.ndns,p.ncwo,p.MSpin)
    elseif p.ipnof==8
        cj12,ck12 = CJCKD8(n,p.no1,p.ndoc,p.nsoc,p.nbeta,p.nalpha,p.ndns,p.ncwo,p.MSpin,p.h_cut)
    end
    cj12[diagind(cj12)] .= 0
    ck12[diagind(ck12)] .= 0            

    return cj12,ck12
end

function der_PNOFi_selector(n,dn_dgamma,p)
    if p.ipnof==5
        Dcj12r,Dck12r = der_CJCKD5(n,dn_dgamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo)
    elseif p.ipnof==7
        Dcj12r,Dck12r = der_CJCKD7(n,p.ista,dn_dgamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo)
    elseif p.ipnof==8
        Dcj12r,Dck12r = der_CJCKD8(n,dn_dgamma,p.no1,p.ndoc,p.nalpha,p.nbeta,p.nv,p.nbf5,p.ndns,p.ncwo,p.MSpin,p.nsoc,p.h_cut)
    end
    for i in 1:p.nv
        for j in 1:p.nbf5
            Dcj12r[j,j,i] = 0
            Dck12r[j,j,i] = 0
        end
    end

    return Dcj12r,Dck12r
end

function CJCKD5(n,no1,ndoc,nsoc,nbeta,nalpha,ndns,ncwo,MSpin)

    # Interpair Electron correlation #

    @tullio cj12[i,j] := 2*n[i]*n[j]
    @tullio ck12[i,j] := n[i]*n[j]

    # Intrapair Electron Correlation #

    if MSpin==0 && nsoc>1
        n_beta_alpha = view(n,nbeta+1:nalpha)
        ck12_beta_alpha = view(ck12,nbeta+1:nalpha,nbeta+1:nalpha)
        @tullio ck12_beta_alpha[i,j] = 2*n_beta_alpha[i]*n_beta_alpha[j]
    end

    for l in 1:ndoc
        ldx = no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + ncwo*(ndoc - l) + 1
        ul = no1 + ndns + ncwo*(ndoc - l + 1)

        cj12[ldx,ll:ul] .= 0
        cj12[ll:ul,ldx] .= 0

        cj12[ll:ul,ll:ul] .= 0

        ck12[ldx,ll:ul] .= sqrt.(n[ldx]*n[ll:ul])
        ck12[ll:ul,ldx] .= sqrt.(n[ldx]*n[ll:ul])

        ck12_ww = view(ck12,ll:ul,ll:ul)
        n_ww = view(n,ll:ul)
        @tullio ck12_ww[i,j] = -sqrt(n_ww[i]*n_ww[j])
    end

    return cj12,ck12

end

function der_CJCKD5(n,dn_dgamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo)

    # Interpair Electron correlation #

    @tullio Dcj12r[i,j,k] := 2*dn_dgamma[i,k]*n[j]
    @tullio Dck12r[i,j,k] := dn_dgamma[i,k]*n[j]

    # Intrapair Electron Correlation

    for l in 1:ndoc
        ldx = no1 + l

        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + ncwo*(ndoc - l) + 1
        ul = no1 + ndns + ncwo*(ndoc - l + 1)

        Dcj12r[ldx,ll:ul,1:nv] .= 0
        Dcj12r[ll:ul,ldx,1:nv] .= 0

        Dcj12r[ll:ul,ll:ul,1:nv] .= 0

        a = max(n[ldx],10^-15)
        b = n[ll:ul]
        b[b.<10^-15] .= 10^-15

        Dck12r_occ_cwo = view(Dck12r,ldx,ll:ul,1:nv)
        Dck12r_cwo_occ = view(Dck12r,ll:ul,ldx,1:nv)
        Dck12r_cwo_cwo = view(Dck12r,ll:ul,ll:ul,1:nv)
        n_cwo = view(n,ll:ul)
        n_occ = n[ldx]
        dn_dgamma_occ = view(dn_dgamma,ldx,1:nv)
        dn_dgamma_cwo = view(dn_dgamma,ll:ul,1:nv)
        @tullio Dck12r_occ_cwo[i,j] = 1/2 * 1/sqrt(a) * dn_dgamma_occ[j] * sqrt(n_cwo[i])
        @tullio Dck12r_cwo_occ[i,j] = 1/2 * 1/sqrt(b[i]) * dn_dgamma_cwo[i,j] * sqrt(n_occ)
        @tullio Dck12r_cwo_cwo[i,j,k] = - 1/2 * 1/sqrt(b[i]) * dn_dgamma_cwo[i,k] * sqrt(n_cwo[j])

    end
    return Dcj12r,Dck12r

end

function CJCKD7(n,ista,no1,ndoc,nsoc,nbeta,nalpha,ndns,ncwo,MSpin)

    if ista==0
        fi = n.*(1 .-n)
        fi[fi.<=0] .= 0
        fi = sqrt.(fi)
    else
        fi = 2*n.*(1 .-n)
    end

    # Interpair Electron correlation #

    @tullio cj12[i,j] := 2*n[i]*n[j]
    @tullio ck12[i,j] := n[i]*n[j] + fi[i]*fi[j]

    # Intrapair Electron Correlation

    if MSpin==0 && nsoc>1
	n_beta_alpha = view(n,nbeta+1:nalpha)
	ck12_beta_alpha = view(ck12,nbeta+1:nalpha,nbeta+1:nalpha)
	@tullio ck12_beta_alpha[i,j] = 2*n_beta_alpha[i]*n_beta_alpha[j]
    end

    for l in 1:ndoc
        ldx = no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + ncwo*(ndoc - l) + 1
        ul = no1 + ndns + ncwo*(ndoc - l + 1)

        cj12[ldx,ll:ul] .= 0
        cj12[ll:ul,ldx] .= 0

        cj12[ll:ul,ll:ul] .= 0

        ck12[ldx,ll:ul] .= sqrt.(n[ldx]*n[ll:ul])
        ck12[ll:ul,ldx] .= sqrt.(n[ldx]*n[ll:ul])

	ck12_ww = view(ck12,ll:ul,ll:ul)
	n_ww = view(n,ll:ul)
	@tullio ck12_ww[i,j] = -sqrt(n_ww[i]*n_ww[j])
    end

    return cj12,ck12

end

function der_CJCKD7(n,ista,dn_dgamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo)

    if ista==0
        fi = n.*(1 .-n)
        fi[fi.<=0] .= 0
        fi = sqrt.(fi)
    else
        fi = 2*n.*(1 .-n)
    end

    dfi_dgamma = zeros(nbf5,nv)
    for i in no1+1:nbf5
        a = max(fi[i],10^-15)
        for k in 1:nv
            if ista==0
		    dfi_dgamma[i,k] = 1/(2*a)*(1-2*n[i])*dn_dgamma[i,k]
            else
		    dfi_dgamma[i,k] = 2*(1-2*n[i])*dn_dgamma[i,k]
            end
	end
    end

    # Interpair Electron correlation #

    @tullio Dcj12r[i,j,k] := 2*dn_dgamma[i,k]*n[j]
    @tullio Dck12r[i,j,k] := dn_dgamma[i,k]*n[j] + dfi_dgamma[i,k]*fi[j]

    # Intrapair Electron Correlation

    for l in 1:ndoc
        ldx = no1 + l

        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
	ll = no1 + ndns + ncwo*(ndoc - l) + 1
        ul = no1 + ndns + ncwo*(ndoc - l + 1)

        Dcj12r[ldx,ll:ul,1:nv] .= 0
        Dcj12r[ll:ul,ldx,1:nv] .= 0

        Dcj12r[ll:ul,ll:ul,1:nv] .= 0

        a = max(n[ldx],10^-15)
        b = n[ll:ul]
        b[b.<10^-15] .= 10^-15

	Dck12r_occ_cwo = view(Dck12r,ldx,ll:ul,1:nv)
	Dck12r_cwo_occ = view(Dck12r,ll:ul,ldx,1:nv)
	Dck12r_cwo_cwo = view(Dck12r,ll:ul,ll:ul,1:nv)
	n_cwo = view(n,ll:ul)
	n_occ = n[ldx]
	dn_dgamma_occ = view(dn_dgamma,ldx,1:nv)
	dn_dgamma_cwo = view(dn_dgamma,ll:ul,1:nv)
	@tullio Dck12r_occ_cwo[i,j] = 1/2 * 1/sqrt(a) * dn_dgamma_occ[j] * sqrt(n_cwo[i])
	@tullio Dck12r_cwo_occ[i,j] = 1/2 * 1/sqrt(b[i]) * dn_dgamma_cwo[i,j] * sqrt(n_occ)
	@tullio Dck12r_cwo_cwo[i,j,k] = - 1/2 * 1/sqrt(b[i]) * dn_dgamma_cwo[i,k] * sqrt(n_cwo[j])

    end
    return Dcj12r,Dck12r

end

function CJCKD8(n,no1,ndoc,nsoc,nbeta,nalpha,ndns,ncwo,MSpin,h_cut)

    nbf5 = size(n)[1]

    #h_cut = p.h_cut#0.02*sqrt(2.0)
    n_d = zeros(nbf5)

    for i in 1:ndoc
        idx = no1 + i
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
	ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        h = 1.0 - n[idx]
        coc = h / h_cut
        arg = -coc^2
        F = exp(arg)  # ! Hd/Hole
        n_d[idx] = n[idx] * F
        n_d[ll:ul] = n[ll:ul] * F  # ROd = RO*Hd/Hole

    end

    n_d12 = sqrt.(n_d)
    fi = n .*(1 .- n)
    fi[fi.<=0] .= 0
    fi = sqrt.(fi)

    # Interpair Electron correlation #

    @tullio cj12[i,j] := 2*n[i]*n[j]
    @tullio ck12[i,j] := n[i]*n[j]

    ck12_beta_cwo = view(ck12,no1+1:nbeta,nalpha+1:nbf5)
    ck12_cwo_beta = view(ck12,nalpha+1:nbf5,no1+1:nbeta)
    ck12_cwo_cwo = view(ck12,nalpha+1:nbf5,nalpha+1:nbf5)
    fi_beta = view(fi,no1+1:nbeta)
    fi_cwo = view(fi,nalpha+1:nbf5)
    @tullio ck12_beta_cwo[i,j] += fi_beta[i]*fi_cwo[j]
    @tullio ck12_cwo_beta[i,j] += fi_cwo[i]*fi_beta[j]
    @tullio ck12_cwo_cwo[i,j] += fi_cwo[i]*fi_cwo[j]

    # Intrapair Electron Correlation

    if(MSpin==0 && nsoc>0)
        ck12[no1+1:nbeta,nbeta+1:nalpha] .+= 0.5*fi[no1+1:nbeta]*0.5
        ck12[nbeta+1:nalpha,no1+1:nbeta] .+= 0.5*0.5*fi[no1+1:nbeta]'
        ck12[nbeta+1:nalpha,nalpha+1:nbf5] .+= 0.5*fi[nalpha+1:nbf5]'
        ck12[nalpha+1:nbf5,nbeta+1:nalpha] .+= fi[nalpha+1:nbf5]*0.5
    end

    if(MSpin==0 && nsoc>1) #then
        ck12[nbeta+1:nalpha,nbeta+1:nalpha] .= 0.5
    end

    n_d12_beta = view(n_d12,no1+1:nbeta)
    n_d12_cwo = view(n_d12,nalpha+1:nbf5)
    n_d_beta = view(n_d,no1+1:nbeta)
    n_d_cwo = view(n_d,nalpha+1:nbf5)
    @tullio ck12_beta_cwo[i,j] += n_d12_beta[i]*n_d12_cwo[j] - n_d_beta[i]*n_d_cwo[j]
    @tullio ck12_cwo_beta[i,j] += n_d12_cwo[i]*n_d12_beta[j] - n_d_cwo[i]*n_d_beta[j]
    @tullio ck12_cwo_cwo[i,j] += -n_d12_cwo[i]*n_d12_cwo[j] - n_d_cwo[i]*n_d_cwo[j]

    for l in 1:ndoc
        ldx = no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + ncwo*(ndoc - l) + 1
        ul = no1 + ndns + ncwo*(ndoc - l + 1)

        cj12[ldx,ll:ul] .= 0
        cj12[ll:ul,ldx] .= 0

        cj12[ll:ul,ll:ul] .= 0

        ck12[ldx,ll:ul] .= sqrt.(n[ldx]*n[ll:ul])
        ck12[ll:ul,ldx] .= sqrt.(n[ldx]*n[ll:ul])

        ck12_ww = view(ck12,ll:ul,ll:ul)
        n_ww = view(n,ll:ul)
        @tullio ck12_ww[i,j] = -sqrt(n_ww[i]*n_ww[j])
    end

    return cj12,ck12

end

function der_CJCKD8(n,dn_dgamma,no1,ndoc,nalpha,nbeta,nv,nbf5,ndns,ncwo,MSpin,nsoc,h_cut)

    nbf5 = size(n)[1]

    #h_cut = p.h_cut#0.02*sqrt(2.0)
    n_d = zeros(nbf5)
    dn_d_dgamma = zeros(nbf5,nv)
    dn_d12_dgamma = zeros(nbf5,nv)

    for i in 1:ndoc
        idx = no1 + i
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        h_idx = 1.0-n[idx]
        coc = h_idx/h_cut
        arg = -coc^2
        F_idx = exp(arg)                # Hd/Hole
        n_d[idx] = n[idx] * F_idx
        n_d[ll:ul] = n[ll:ul] * F_idx      # n_d = RO*Hd/Hole
        dn_d_dgamma[idx,1:nv] .= F_idx*dn_dgamma[idx,1:nv] .* (1-n[idx]*(- 2*coc/h_cut))
        dn_d_dgamma[ll:ul,1:nv] .= F_idx*dn_dgamma[ll:ul,1:nv]
	dn_d_dgamma_cwo_v = view(dn_d_dgamma,ll:ul,1:nv)
	dn_dgamma_v = view(dn_dgamma,idx,1:nv)
	n_cwo = view(n,ll:ul)
	@tullio dn_d_dgamma_cwo_v[i,j] += F_idx*(2*coc/h_cut) * n_cwo[i] * dn_dgamma_v[j]
    end

    n_d12 = sqrt.(n_d)
    dn_d12_dgamma = 0.5 * dn_d_dgamma ./ n_d12

    fi = n .*(1 .- n)
    fi[fi.<=0] .= 0
    fi = sqrt.(fi)

    dfi_dgamma = zeros(nbf5,nv)
    for i in no1+1:nbf5
        a = max(fi[i],10^-15)
        for k in 1:nv
            dfi_dgamma[i,k] = 1/(2*a)*(1-2*n[i])*dn_dgamma[i,k]
        end
    end


    # Interpair Electron correlation #

    @tullio Dcj12r[i,j,k] := 2*dn_dgamma[i,k]*n[j]
    @tullio Dck12r[i,j,k] := dn_dgamma[i,k]*n[j]

    Dck12r_beta_cwo = view(Dck12r,no1+1:nbeta,nalpha+1:nbf5,1:nv)
    Dck12r_cwo_beta = view(Dck12r,nalpha+1:nbf5,no1+1:nbeta,1:nv)
    Dck12r_cwo_cwo = view(Dck12r,nalpha+1:nbf5,nalpha+1:nbf5,1:nv)
    fi_beta = view(fi,no1+1:nbeta)
    fi_cwo = view(fi,nalpha+1:nbf5)
    dfi_beta = view(dfi_dgamma,no1+1:nbeta,1:nv)
    dfi_cwo = view(dfi_dgamma,nalpha+1:nbf5,1:nv)
    @tullio Dck12r_beta_cwo[i,j,k] += dfi_beta[i,k]*fi_cwo[j]
    @tullio Dck12r_cwo_beta[i,j,k] += fi_beta[j]*dfi_cwo[i,k]
    @tullio Dck12r_cwo_cwo[i,j,k] += fi_cwo[j]*dfi_cwo[i,k]

    if(MSpin==0 && nsoc>0)
        Dck12r_beta_alpha = view(Dck12r,no1+1:nbeta,nbeta+1:nalpha,1:nv)
        Dck12r_cwo_alpha = view(Dck12r,nalpha+1:nbf5,nbeta+1:nalpha,1:nv)
        dfi_cwo = view(dfi_dgamma,nalpha+1:nbf5,1:nv)
	@tullio Dck12r_beta_alpha[i,j,k] += 0.5*dfi_beta[i,k]*0.5
	@tullio Dck12r_cwo_alpha[i,j,k] += dfi_cwo[i,k]*0.5
    end

    if(MSpin==0 && nsoc>1)
        Dck12r[nbeta+1:nalpha, nbeta+1:nalpha, 1:nv] .= 0.0
    end

    n_d12_beta = view(n_d12,no1+1:nbeta)
    n_d12_cwo = view(n_d12,nalpha+1:nbf5)
    n_d_beta = view(n_d,no1+1:nbeta)
    n_d_cwo = view(n_d,nalpha+1:nbf5)
    dn_d12_beta = view(dn_d12_dgamma,no1+1:nbeta,1:nv)
    dn_d12_cwo = view(dn_d12_dgamma,nalpha+1:nbf5,1:nv)
    dn_d_beta = view(dn_d_dgamma,no1+1:nbeta,1:nv)
    dn_d_cwo = view(dn_d_dgamma,nalpha+1:nbf5,1:nv)
    @tullio Dck12r_beta_cwo[i,j,k] += dn_d12_beta[i,k]*n_d12_cwo[j] - dn_d_beta[i,k]*n_d_cwo[j]
    @tullio Dck12r_cwo_beta[i,j,k] += n_d12_beta[j]*dn_d12_cwo[i,k] - n_d_beta[j]*dn_d_cwo[i,k]
    @tullio Dck12r_cwo_cwo[i,j,k] += -n_d12_cwo[j]*dn_d12_cwo[i,k] - n_d_cwo[j]*dn_d_cwo[i,k]

    # Intrapair Electron Correlation

    for l in 1:ndoc
        ldx = no1 + l

        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = no1 + ndns + ncwo*(ndoc - l) + 1
        ul = no1 + ndns + ncwo*(ndoc - l + 1)

        Dcj12r[ldx,ll:ul,1:nv] .= 0
        Dcj12r[ll:ul,ldx,1:nv] .= 0

        Dcj12r[ll:ul,ll:ul,1:nv] .= 0

        a = max(n[ldx],10^-15)
        b = n[ll:ul]
        b[b.<10^-15] .= 10^-15

        Dck12r_occ_cwo = view(Dck12r,ldx,ll:ul,1:nv)
        Dck12r_cwo_occ = view(Dck12r,ll:ul,ldx,1:nv)
        Dck12r_cwo_cwo = view(Dck12r,ll:ul,ll:ul,1:nv)
        n_cwo = view(n,ll:ul)
        n_occ = n[ldx]
        dn_dgamma_occ = view(dn_dgamma,ldx,1:nv)
        dn_dgamma_cwo = view(dn_dgamma,ll:ul,1:nv)
        @tullio Dck12r_occ_cwo[i,j] = 1/2 * 1/sqrt(a) * dn_dgamma_occ[j] * sqrt(n_cwo[i])
        @tullio Dck12r_cwo_occ[i,j] = 1/2 * 1/sqrt(b[i]) * dn_dgamma_cwo[i,j] * sqrt(n_occ)
        @tullio Dck12r_cwo_cwo[i,j,k] = - 1/2 * 1/sqrt(b[i]) * dn_dgamma_cwo[i,k] * sqrt(n_cwo[j])

    end

    return Dcj12r,Dck12r

end

function ocupacion(gamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo,HighSpin,occ_method)
    if occ_method == "Trigonometric"
        n,dn_dgamma = ocupacion_trigonometric(gamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo,HighSpin)
    elseif occ_method == "Softmax"
        n,dn_dgamma = ocupacion_softmax(gamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo,HighSpin)
    end
    return n,dn_dgamma
end

function ocupacion_trigonometric(gamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo,HighSpin)

    n = zeros(nbf5)
    dni_dgammai = zeros(nbf5)
    dn_dgamma = zeros(nbf5,nv)

    n[1:no1] .= 1                                              # [1,no1]

    n[no1+1:no1+ndoc] .= 1/2 * (1 .+ (cos.(gamma[1:ndoc])).^2)     # (no1,no1+ndoc]
    dni_dgammai[no1+1:no1+ndoc] .= - 1/2 * sin.(2*gamma[1:ndoc])

    if !HighSpin
        n[no1+ndoc+1:no1+ndns] .= 0.5   # (no1+ndoc,no1+ndns]
    elseif HighSpin
        n[no1+ndoc+1:no1+ndns] .= 1.0   # (no1+ndoc,no1+ndns]
    end

    dn_dgamma = zeros(nbf5,nv)
    h = 1 .- n

    for i in 1:ndoc
        ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        n[ll:ul] .= h[no1+i]
        for iw in 1:ncwo-1
            n[ll+iw-1] *= sin(gamma[ndoc+(ncwo-1)*(i-1)+iw])^2
            n[ll+iw:ul] .*= cos(gamma[ndoc+(ncwo-1)*(i-1)+iw])^2
        end
    end

    for i in 1:ndoc
        # dn_g/dgamma_g
        dn_dgamma[no1+i,i] = dni_dgammai[no1+i]

        # dn_pi/dgamma_g
        ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        dn_dgamma[ll:ul,i] .= -dni_dgammai[no1+i]
        for iw in 1:ncwo-1
            dn_dgamma[ll+iw-1,i] *= sin(gamma[ndoc+(ncwo-1)*(i-1)+iw])^2
    	dn_dgamma[ll+iw:ul,i] .*= cos(gamma[ndoc+(ncwo-1)*(i-1)+iw])^2
        end

        # dn_pi/dgamma_pj (j<i)
        for iw in 1:ncwo-1
    	dn_dgamma[ll+iw:ul,ndoc+(ncwo-1)*(i-1)+iw] .= n[no1+i] - 1
            for ip in ll+iw:ul
                for jw in 1:ip-ll
                    if jw==iw
                        dn_dgamma[ip,ndoc+(ncwo-1)*(i-1)+iw] *= sin(2*gamma[ndoc+(ncwo-1)*(i-1)+jw])
                    else
    		    dn_dgamma[ip,ndoc+(ncwo-1)*(i-1)+iw] *= cos(gamma[ndoc+(ncwo-1)*(i-1)+jw])^2
                    end
                end
                if ip-ll+1<=ncwo-1
                    dn_dgamma[ip,ndoc+(ncwo-1)*(i-1)+iw] *= sin(gamma[ndoc+(ncwo-1)*(i-1)+(ip-ll+1)])^2
                end
            end
        end

        # dn_pi/dgamma_i
        for iw in 1:ncwo-1
            dn_dgamma[ll+iw-1,ndoc+(ncwo-1)*(i-1)+iw] = 1 - n[no1+i]
            for jw in 1:iw
                if jw==iw
                    dn_dgamma[ll+iw-1,ndoc+(ncwo-1)*(i-1)+iw] *= sin(2*gamma[ndoc+(ncwo-1)*(i-1)+jw])
                else
                    dn_dgamma[ll+iw-1,ndoc+(ncwo-1)*(i-1)+iw] *= cos(gamma[ndoc+(ncwo-1)*(i-1)+jw])^2
    	    end
    	end
        end
    end

   return n,dn_dgamma
end

function ocupacion_softmax(gamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo,HighSpin)

    n = zeros(nbf5)
    dn_dgamma = zeros(nbf5,nv)

    n[1:no1] .= 1                                              # [1,no1]
    if !HighSpin
        n[no1+ndoc+1:no1+ndns] .= 0.5   # (no1+ndoc,no1+ndns]
    elseif HighSpin
        n[no1+ndoc+1:no1+ndns] .= 1.0   # (no1+ndoc,no1+ndns]
    end

    exp_gamma = exp.(gamma)

    for i in 1:ndoc

        ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        llg = ndoc + ncwo*(ndoc - i) + 1
        ulg = ndoc + ncwo*(ndoc - i + 1)

        sum_exp = exp_gamma[i] + sum(exp_gamma[llg:ulg])
	
        n[i] = exp_gamma[i] / sum_exp 
        n[ll:ul] .= exp_gamma[llg:ulg] ./ sum_exp 

        for j in ll:ul
            for k in llg:ulg
                jg = j - no1 - ndns + ndoc
                dn_dgamma[j,k] = -exp_gamma[jg]*exp_gamma[k]/sum_exp^2
            end
        end
        for jg in llg:ulg
            dn_dgamma[i,jg] = -exp_gamma[i]*exp_gamma[jg]/sum_exp^2
        end
        for j in ll:ul
            jg = j - no1 - ndns + ndoc
            dn_dgamma[j,i] = -exp_gamma[i]*exp_gamma[jg]/sum_exp^2
        end

        dn_dgamma[i,i] = exp_gamma[i]*(sum_exp - exp_gamma[i])/sum_exp ^2
        for j in ll:ul
            jg = j - no1 - ndns + ndoc
            dn_dgamma[j,jg] = exp_gamma[jg]*(sum_exp - exp_gamma[jg])/sum_exp ^2
        end
    end
    return n,dn_dgamma
end

function calce(n,cj12,ck12,J_MO,K_MO,H_core,p)

    E = 0

    n_beta = view(n,1:p.nbeta)
    n_alpha = view(n,p.nbeta+1:p.nalpha)
    n_nbf5 = view(n,p.nalpha+1:p.nbf5)
    H_beta = view(H_core,1:p.nbeta)
    H_alpha= view(H_core,p.nbeta+1:p.nalpha)
    H_nbf5 = view(H_core,p.nalpha+1:p.nbf5)
    J_MO_beta = view(J_MO,1:p.nbeta,1:p.nbeta)
    J_MO_nbf5 = view(J_MO,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
    if p.MSpin==0

        # 2H + J
        @tullio E += n_beta[i]*(2*H_beta[i] + J_MO_beta[i,i])
	if(p.nsoc>0)
            @tullio E += n_alpha[i]*2*H_alpha[i]
        end
	@tullio E += n_nbf5[i]*(2*H_nbf5[i]+J_MO_nbf5[i,i])

        #C^J JMO
	@tullio E += cj12[i,j]*J_MO[j,i]

        #C^K KMO
	@tullio E += -ck12[i,j]*K_MO[j,i]

    elseif p.MSpin!=0

        # 2H + J
        @tullio E += n_beta[i]*(2*H_beta[i] + J_MO_beta[i,i])
	@tullio E += n_alpha[i]*2*H_alpha[i]
	@tullio E += n_nbf5[i]*(2*H_nbf5[i]+J_MO_nbf5[i,i])

        #C^J JMO
        cj12_beta_beta = view(cj12,1:p.nbeta,1:p.nbeta)
        cj12_beta_nbf5 = view(cj12,1:p.nbeta,p.nalpha+1:p.nbf5)
        cj12_nbf5_beta = view(cj12,p.nalpha+1:p.nbf5,1:p.nbeta)
        cj12_nbf5_nbf5 = view(cj12,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
        J_MO_beta_beta = view(J_MO,1:p.nbeta,1:p.nbeta)
        J_MO_nbf5_beta = view(J_MO,p.nalpha+1:p.nbf5,1:p.nbeta)
        J_MO_beta_nbf5 = view(J_MO,1:p.nbeta,p.nalpha+1:p.nbf5)
        J_MO_nbf5_nbf5 = view(J_MO,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
        J_MO_alpha_beta = view(J_MO,p.nbeta+1:p.nalpha,1:p.nbeta)
        J_MO_alpha_alpha = view(J_MO,p.nbeta+1:p.nalpha,p.nbeta+1:p.nalpha)
        J_MO_alpha_nbf5 = view(J_MO,p.nbeta+1:p.nalpha,p.nalpha+1:p.nbf5)
	@tullio E += cj12_beta_beta[i,j]*J_MO_beta_beta[j,i]
	@tullio E += cj12_beta_nbf5[i,j]*J_MO_nbf5_beta[j,i]
	@tullio E += cj12_nbf5_beta[i,j]*J_MO_beta_nbf5[j,i]
	@tullio E += cj12_nbf5_nbf5[i,j]*J_MO_nbf5_nbf5[j,i]

        #C^K KMO
        ck12_beta_beta = view(ck12,1:p.nbeta,1:p.nbeta)
        ck12_beta_nbf5 = view(ck12,1:p.nbeta,p.nalpha+1:p.nbf5)
        ck12_nbf5_beta = view(ck12,p.nalpha+1:p.nbf5,1:p.nbeta)
        ck12_nbf5_nbf5 = view(ck12,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
        K_MO_beta_beta = view(K_MO,1:p.nbeta,1:p.nbeta)
        K_MO_nbf5_beta = view(K_MO,p.nalpha+1:p.nbf5,1:p.nbeta)
        K_MO_beta_nbf5 = view(K_MO,1:p.nbeta,p.nalpha+1:p.nbf5)
        K_MO_nbf5_nbf5 = view(K_MO,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
        K_MO_alpha_beta = view(K_MO,p.nbeta+1:p.nalpha,1:p.nbeta)
        K_MO_alpha_alpha = view(K_MO,p.nbeta+1:p.nalpha,p.nbeta+1:p.nalpha)
	K_MO_alpha_nbf5 = view(K_MO,p.nbeta+1:p.nalpha,p.nalpha+1:p.nbf5)
	@tullio E += -ck12_beta_beta[i,j]*K_MO_beta_beta[j,i]
	@tullio E += -ck12_beta_nbf5[i,j]*K_MO_nbf5_beta[j,i]
	@tullio E += -ck12_nbf5_beta[i,j]*K_MO_beta_nbf5[j,i]
	@tullio E += -ck12_nbf5_nbf5[i,j]*K_MO_nbf5_nbf5[j,i]

        #n JMO
	@tullio E += 2*n_beta[i]*J_MO_alpha_beta[j,i]
	@tullio E += 2*n_nbf5[i]*J_MO_alpha_nbf5[j,i]
	@tullio E += 0.5*n_alpha[i]*J_MO_alpha_alpha[j,i]
	@tullio E += -0.5*n_alpha[i]*J_MO_alpha_alpha[i,i]

        #n KMO
	@tullio E += -n_beta[i]*K_MO_alpha_beta[j,i]
	@tullio E += -n_nbf5[i]*K_MO_alpha_nbf5[j,i]
	@tullio E += -n_alpha[i]*K_MO_alpha_alpha[j,i]
	@tullio E += n_alpha[i]*K_MO_alpha_alpha[i,i]
    end

    return E
end

function calcocce(gamma,J_MO,K_MO,H_core,p)

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
    cj12,ck12 = PNOFi_selector(n,p)

    E = calce(n,cj12,ck12,J_MO,K_MO,H_core,p)

    return E

end

function calcoccg(gamma,J_MO,K_MO,H_core,p)

    grad = zeros(p.nv)
    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
    Dcj12r,Dck12r = der_PNOFi_selector(n,dn_dgamma,p)

    if p.MSpin==0

        # dn_dgamma (2H+J)
	dn_dgamma_beta = view(dn_dgamma,p.no1+1:p.nbeta,1:p.nv)
	dn_dgamma_nbf5 = view(dn_dgamma,p.nalpha+1:p.nbf5,1:p.nv)
	H_core_beta = view(H_core,p.no1+1:p.nbeta)
	H_core_nbf5 = view(H_core,p.nalpha+1:p.nbf5)
	J_MO_beta = view(J_MO,p.no1+1:p.nbeta,p.no1+1:p.nbeta)
	J_MO_nbf5 = view(J_MO,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
	@tullio grad[k] += dn_dgamma_beta[i,k]*2*H_core_beta[i]# + J_MO_beta[i,i])
	@tullio grad[k] += dn_dgamma_beta[i,k]*J_MO_beta[i,i]
        #grad += np.einsum('ik,i->k',dn_dgamma[p.no1:p.nbeta,:p.nv],2*H_core[p.no1:p.nbeta]+np.diagonal(J_MO)[p.no1:p.nbeta],optimize=True) # [0,Nbeta]
	@tullio grad[k] += dn_dgamma_nbf5[i,k]*2*H_core_nbf5[i]# + J_MO_nbf5[i,i])
	@tullio grad[k] += dn_dgamma_nbf5[i,k]*J_MO_nbf5[i,i]

        # 2 dCJ_dgamma J_MO
	Dcj12r_beta = view(Dcj12r,p.no1+1:p.nbeta,1:p.nbf5,1:p.nv)
	Dcj12r_nbf5 = view(Dcj12r,p.nalpha+1:p.nbf5,1:p.nbf5,1:p.nv)
	J_MO_beta = view(J_MO,1:p.nbf5,p.no1+1:p.nbeta)
	J_MO_nbf5 = view(J_MO,1:p.nbf5,p.nalpha+1:p.nbf5)
	@tullio grad[k] += 2*Dcj12r_beta[i,j,k]*J_MO_beta[j,i]
        #grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.no1:p.nbeta,:p.nbf5,:p.nv],J_MO[:p.nbf5,p.no1:p.nbeta],optimize=True)

	@tullio grad[k] += 2*Dcj12r_nbf5[i,j,k]*J_MO_nbf5[j,i]
        #grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.nalpha:p.nbf5,:p.nbf5,:p.nv],J_MO[:p.nbf5,p.nalpha:p.nbf5],optimize=True)

        # -2 dCK_dgamma K_MO
	Dck12r_beta = view(Dck12r,p.no1+1:p.nbeta,1:p.nbf5,1:p.nv)
	Dck12r_nbf5 = view(Dck12r,p.nalpha+1:p.nbf5,1:p.nbf5,1:p.nv)
	K_MO_beta = view(K_MO,1:p.nbf5,p.no1+1:p.nbeta)
	K_MO_nbf5 = view(K_MO,1:p.nbf5,p.nalpha+1:p.nbf5)
	@tullio grad[k] += -2*Dck12r_beta[i,j,k]*K_MO_beta[j,i]
        #grad -= 2*np.einsum('ijk,ji->k',Dck12r[p.no1:p.nbeta,:p.nbf5,:p.nv],K_MO[:p.nbf5,p.no1:p.nbeta],optimize=True)

	@tullio grad[k] += -2*Dck12r_nbf5[i,j,k]*K_MO_nbf5[j,i]
        #grad -= 2*np.einsum('ijk,ji->k',Dck12r[p.nalpha:p.nbf5,:p.nbf5,:p.nv],K_MO[:p.nbf5,p.nalpha:p.nbf5],optimize=True)

    elseif p.MSpin!=0
#
#        # dn_dgamma (2H+J)
        dn_dgamma_beta = view(dn_dgamma,p.no1+1:p.nbeta,1:p.nv)
        dn_dgamma_nbf5 = view(dn_dgamma,p.nalpha+1:p.nbf5,1:p.nv)
        H_core_beta = view(H_core,p.no1+1:p.nbeta)
        H_core_nbf5 = view(H_core,p.nalpha+1:p.nbf5)
        J_MO_beta = view(J_MO,p.no1+1:p.nbeta,p.no1+1:p.nbeta)
        J_MO_nbf5 = view(J_MO,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
        @tullio grad[k] += dn_dgamma_beta[i,k]*2*H_core_beta[i]# + J_MO_beta[i,i])
        @tullio grad[k] += dn_dgamma_beta[i,k]*J_MO_beta[i,i]
        #grad += np.einsum('ik,i->k',dn_dgamma[p.no1:p.nbeta,:p.nv],2*H_core[p.no1:p.nbeta]+np.diagonal(J_MO)[p.no1:p.nbeta],optimize=True) # [0,Nbeta]
        @tullio grad[k] += dn_dgamma_nbf5[i,k]*2*H_core_nbf5[i]# + J_MO_nbf5[i,i])
        @tullio grad[k] += dn_dgamma_nbf5[i,k]*J_MO_nbf5[i,i]
#        grad += np.einsum('ik,i->k',dn_dgamma[p.no1:p.nbeta,:p.nv],2*H_core[p.no1:p.nbeta]+np.diagonal(J_MO)[p.no1:p.nbeta],optimize=True) # [0,Nbeta]
#        grad += np.einsum('ik,i->k',dn_dgamma[p.nalpha:p.nbf5,:p.nv],2*H_core[p.nalpha:p.nbf5]+np.diagonal(J_MO)[p.nalpha:p.nbf5],optimize=True) # [Nalpha,Nbf5]
#
#        # 2 dCJ_dgamma J_MO
        Dcj12r_beta_beta = view(Dcj12r,p.no1+1:p.nbeta,1:p.nbeta,1:p.nv)
        Dcj12r_nbf5_beta = view(Dcj12r,p.nalpha+1:p.nbf5,1:p.nbeta,1:p.nv)
        Dcj12r_beta_nbf5 = view(Dcj12r,p.no1+1:p.nbeta,p.nalpha+1:p.nbf5,1:p.nv)
        Dcj12r_nbf5_nbf5 = view(Dcj12r,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5,1:p.nv)
        J_MO_beta_beta = view(J_MO,1:p.nbeta,p.no1+1:p.nbeta)
        J_MO_beta_nbf5 = view(J_MO,1:p.nbeta,p.nalpha+1:p.nbf5)
        J_MO_nbf5_beta = view(J_MO,p.nalpha+1:p.nbf5,p.no1+1:p.nbeta)
        J_MO_nbf5_nbf5 = view(J_MO,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
        @tullio grad[k] += 2*Dcj12r_beta_beta[i,j,k]*J_MO_beta_beta[j,i]
        @tullio grad[k] += 2*Dcj12r_beta_nbf5[i,j,k]*J_MO_nbf5_beta[j,i]
        @tullio grad[k] += 2*Dcj12r_nbf5_beta[i,j,k]*J_MO_beta_nbf5[j,i]
        @tullio grad[k] += 2*Dcj12r_nbf5_nbf5[i,j,k]*J_MO_nbf5_nbf5[j,i]
#        grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.no1:p.nbeta,:p.nbeta,:p.nv],J_MO[:p.nbeta,p.no1:p.nbeta],optimize=True)
#        grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.no1:p.nbeta,p.nalpha:p.nbf5,:p.nv],J_MO[p.nalpha:p.nbf5,p.no1:p.nbeta],optimize=True)
#        grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.nalpha:p.nbf5,:p.nbeta,:p.nv],J_MO[:p.nbeta,p.nalpha:p.nbf5],optimize=True)
#        grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.nalpha:p.nbf5,p.nalpha:p.nbf5,:p.nv],J_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5],optimize=True)
#
#        # -2 dCK_dgamma K_MO
        Dck12r_beta_beta = view(Dck12r,p.no1+1:p.nbeta,1:p.nbeta,1:p.nv)
        Dck12r_nbf5_beta = view(Dck12r,p.nalpha+1:p.nbf5,1:p.nbeta,1:p.nv)
        Dck12r_beta_nbf5 = view(Dck12r,p.no1+1:p.nbeta,p.nalpha+1:p.nbf5,1:p.nv)
        Dck12r_nbf5_nbf5 = view(Dck12r,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5,1:p.nv)
        K_MO_beta_beta = view(K_MO,1:p.nbeta,p.no1+1:p.nbeta)
        K_MO_beta_nbf5 = view(K_MO,1:p.nbeta,p.nalpha+1:p.nbf5)
        K_MO_nbf5_beta = view(K_MO,p.nalpha+1:p.nbf5,p.no1+1:p.nbeta)
        K_MO_nbf5_nbf5 = view(K_MO,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
        @tullio grad[k] += -2*Dck12r_beta_beta[i,j,k]*K_MO_beta_beta[j,i]
        @tullio grad[k] += -2*Dck12r_beta_nbf5[i,j,k]*K_MO_nbf5_beta[j,i]
        @tullio grad[k] += -2*Dck12r_nbf5_beta[i,j,k]*K_MO_beta_nbf5[j,i]
        @tullio grad[k] += -2*Dck12r_nbf5_nbf5[i,j,k]*K_MO_nbf5_nbf5[j,i]
#        grad -= 2*np.einsum('ijk,ji->k',Dck12r[p.no1:p.nbeta,:p.nbeta,:p.nv],K_MO[:p.nbeta,p.no1:p.nbeta],optimize=True)
#        grad -= 2*np.einsum('ijk,ji->k',Dck12r[p.no1:p.nbeta,p.nalpha:p.nbf5,:p.nv],K_MO[p.nalpha:p.nbf5,p.no1:p.nbeta],optimize=True)
#        grad -= 2*np.einsum('ijk,ji->k',Dck12r[p.nalpha:p.nbf5,:p.nbeta,:p.nv],K_MO[:p.nbeta,p.nalpha:p.nbf5],optimize=True)
#        grad -= 2*np.einsum('ijk,ji->k',Dck12r[p.nalpha:p.nbf5,p.nalpha:p.nbf5,:p.nv],K_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5],optimize=True)
#
#
#        # 2 dn_dgamma J_MO
        dn_dgamma_beta = view(dn_dgamma,p.no1+1:p.nbeta,1:p.nv)
        dn_dgamma_nbf5 = view(dn_dgamma,p.nalpha+1:p.nbf5,1:p.nv)
	J_MO_alpha_beta = view(J_MO,p.nbeta+1:p.nalpha,p.no1+1:p.nbeta)
	J_MO_alpha_nbf5 = view(J_MO,p.nbeta+1:p.nalpha,p.nalpha+1:p.nbf5)
	@tullio grad[k] += 2*dn_dgamma_beta[j,k]*J_MO_alpha_beta[i,j]
	@tullio grad[k] += 2*dn_dgamma_nbf5[j,k]*J_MO_alpha_nbf5[i,j]
#        grad += 2*np.einsum('jk,ij->k',dn_dgamma[p.no1:p.nbeta,:p.nv],J_MO[p.nbeta:p.nalpha,p.no1:p.nbeta],optimize=True)
#        grad += 2*np.einsum('jk,ij->k',dn_dgamma[p.nalpha:p.nbf5,:p.nv],J_MO[p.nbeta:p.nalpha,p.nalpha:p.nbf5],optimize=True)
#
#        # - dn_dgamma K_MO
	K_MO_alpha_beta = view(K_MO,p.nbeta+1:p.nalpha,p.no1+1:p.nbeta)
	K_MO_alpha_nbf5 = view(K_MO,p.nbeta+1:p.nalpha,p.nalpha+1:p.nbf5)
	@tullio grad[k] += -dn_dgamma_beta[j,k]*K_MO_alpha_beta[i,j]
	@tullio grad[k] += -dn_dgamma_nbf5[j,k]*K_MO_alpha_nbf5[i,j]
#        grad -= np.einsum('jk,ij->k',dn_dgamma[p.no1:p.nbeta,:p.nv],K_MO[p.nbeta:p.nalpha,p.no1:p.nbeta],optimize=True)
#        grad -= np.einsum('jk,ij->k',dn_dgamma[p.nalpha:p.nbf5,:p.nv],K_MO[p.nbeta:p.nalpha,p.nalpha:p.nbf5],optimize=True)
    end
    return grad
end

function calcorbe(y,n,cj12,ck12,C,H,I,b_mnl,p)

    Cnew = rotate_orbital(y,C,p)
    J_MO,K_MO,H_core = computeJKH_MO(Cnew,H,I,b_mnl,p)
    E = calce(n,cj12,ck12,J_MO,K_MO,H_core,p)

    return E

end

function calcorbg(y,n,cj12,ck12,C,H,I,b_mnl,pa)

    Cnew = rotate_orbital(y,C,pa)

    elag,Hmat = compute_Lagrange2(Cnew,n,H,I,b_mnl,cj12,ck12,pa)

    grad = 4*elag - 4*elag'
    grads = zeros(pa.nvar)
    n = 1
    for i in 1:pa.nbf5
        for j in i+1:pa.nbf
            grads[n] = grad[i,j]
            n += 1
        end
    end

    return grads

end

function calcorbeg(F,G,y,n,cj12,ck12,C,H,I,b_mnl,pa)
    
    Cnew = rotate_orbital(y,C,pa)

    elag,Hmat = compute_Lagrange2(Cnew,n,H,I,b_mnl,cj12,ck12,pa)

    if G!==nothing
        grad = 4*elag - 4*elag'
        grads = zeros(pa.nvar)
        nn = 1
        for i in 1:pa.nbf5
            for j in i+1:pa.nbf
                grads[nn] = grad[i,j]
                nn += 1
            end
        end
        G .= grads
    end
    if F!==nothing
        E = computeE_elec(Hmat,n,elag,pa)
        return E
    end

end

function calccombe(x,C,H,I,b_mnl,p)

    y = x[1:p.nvar]
    gamma = x[p.nvar+1:end]

    Cnew = rotate_orbital(y,C,p)

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
    cj12,ck12 = PNOFi_selector(n,p)

    J_MO,K_MO,H_core = computeJKH_MO(Cnew,H,I,b_mnl,p)

    E = calce(n,cj12,ck12,J_MO,K_MO,H_core,p)

    return E
end

function calccombg(x,C,H,I,b_mnl,p)

    y = x[1:p.nvar]
    gamma = x[p.nvar+1:end]

    Cnew = rotate_orbital(y,C,p)

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
    cj12,ck12 = PNOFi_selector(n,p)

    J_MO,K_MO,H_core = computeJKH_MO(Cnew,H,I,b_mnl,p)

    grad = zeros(p.nvar + p.nv)

    grad[1:p.nvar] = calcorbg(y,n,cj12,ck12,C,H,I,b_mnl,p)
    grad[p.nvar+1:end] = calcoccg(gamma,J_MO,K_MO,H_core,p)

    return grad
end

function calccombeg(F,G,x,C,H,I,b_mnl,p)

    y = x[1:p.nvar]
    gamma = x[p.nvar+1:end]

    Cnew = rotate_orbital(y,C,p)

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
    cj12,ck12 = PNOFi_selector(n,p)

    J_MO,K_MO,H_core = computeJKH_MO(Cnew,H,I,b_mnl,p)

    if G!==nothing
        grad = zeros(p.nvar + p.nv)

        grad[p.nvar+1:end] = calcoccg(gamma,J_MO,K_MO,H_core,p)
        if F!==nothing
	    gv = view(grad,1:p.nvar)
            E = calcorbeg(F,gv,y,n,cj12,ck12,C,H,I,b_mnl,p)
	    G .= grad
	    return E
        else
            grad[1:p.nvar] = calcorbg(y,n,cj12,ck12,C,H,I,b_mnl,p)
	    G .= grad
	end
    end
    if F!==nothing
        E = calcocce(gamma,J_MO,K_MO,H_core,p)
	return E
    end

end

function compute_2RDM(pp,n)

    # PNOF5
    # Interpair Electron correlation
    Id = 1. * Matrix(I, pp.nbf5, pp.nbf5)

    inter = n .* n'
    intra = zeros(pp.nbf5,pp.nbf5)

    # Intrapair Electron Correlation
    for l in 1:pp.ndoc
        ldx = pp.no1 + l
        # inicio y fin de los orbitales acoplados a los fuertemente ocupados
        ll = pp.no1 + pp.ndns + pp.ncwo*(pp.ndoc - l) + 1
        ul = pp.no1 + pp.ndns + pp.ncwo*(pp.ndoc - l + 1)

        inter[ldx,ldx] = 0
        inter[ldx,ll:ul] .= 0
        inter[ll:ul,ldx] .= 0
        inter[ll:ul,ll:ul] .= 0

        intra[ldx,ldx] = sqrt(n[ldx]*n[ldx])
        intra[ldx,ll:ul] .= -sqrt.(n[ldx]*n[ll:ul])
        intra[ll:ul,ldx] .= -sqrt.(n[ldx]*n[ll:ul])

        intra_ww = view(intra,ll:ul,ll:ul)
        n_ww = view(n,ll:ul)
        @tullio intra_ww[i,j] = sqrt(n_ww[i]*n_ww[j])
    end

    for i in pp.nbeta+1:pp.nalpha
        inter[i,i] = 0
    end

    @tullio Daa[p,q,r,t] :=  inter[p,q]*Id[p,r]*Id[q,t]
    @tullio Daa[p,q,r,t] += -inter[p,q]*Id[p,t]*Id[q,r]
    @tullio Dab[p,q,r,t] :=  inter[p,q]*Id[p,r]*Id[q,t]
    @tullio Dab[p,q,r,t] +=  intra[p,r]*Id[p,q]*Id[r,t]

    # PNOF7
    if(pp.ipnof == 7 || pp.ipnof==8)
        fi = n .* (1 .- n)
        fi[fi .<= 0] .= 0
        fi = sqrt.(fi)
        Pi_s = fi .* fi'
        # Intrapair Electron Correlation
        for l in 1:pp.ndoc
            ldx = pp.no1 + l
            # inicio y fin de los orbitales acoplados a los fuertemente ocupados
            ll = pp.no1 + pp.ndns + pp.ncwo*(pp.ndoc - l) + 1
            ul = pp.no1 + pp.ndns + pp.ncwo*(pp.ndoc - l + 1)

            Pi_s[ldx,ldx] = 0
            Pi_s[ldx,ll:ul] .= 0
            Pi_s[ll:ul,ldx] .= 0
            Pi_s[ll:ul,ll:ul] .= 0
	end
        for i in pp.nbeta+1:pp.nalpha
            Pi_s[i,i] = 0
        end

        inter2 = 0*Pi_s
        inter2[pp.nbeta+1:pp.nalpha,pp.nbeta+1:pp.nalpha] = Pi_s[pp.nbeta+1:pp.nalpha,pp.nbeta+1:pp.nalpha]
        Pi_s[pp.nbeta+1:pp.nalpha,pp.nbeta+1:pp.nalpha] .= 0

	@tullio Dab[p,q,r,t] += -inter2[p,r]*Id[p,t]*Id[q,r]

        if(pp.ipnof==8)
            Pi_s[1:pp.nbeta,1:pp.nbeta] .= 0
            Pi_s[1:pp.nbeta,pp.nbeta+1:pp.nalpha] *= 0.5
            Pi_s[pp.nbeta+1:pp.nalpha,1:pp.nbeta] *= 0.5
        end

	@tullio Dab[p,q,r,t] += -Pi_s[p,r]*Id[p,q]*Id[r,t]

        if(pp.ipnof==8)

            h_cut = 0.02*sqrt(2.0)
	    n_d = zeros(size(n)[1])

            for i in 1:pp.ndoc
                idx = pp.no1 + i
                # inicio y fin de los orbitales acoplados a los fuertemente ocupados
                ll = pp.no1 + pp.ndns + pp.ncwo*(pp.ndoc - i) + 1
                ul = pp.no1 + pp.ndns + pp.ncwo*(pp.ndoc - i + 1)

                h = 1.0 - n[idx]
                coc = h / h_cut
                arg = -coc^2
                F = exp(arg)  # ! Hd/Hole
                n_d[idx] = n[idx] * F
                n_d[ll:ul] = n[ll:ul] * F  # ROd = RO*Hd/Hole
            end

            n_d12 = sqrt.(n_d)

            inter = n_d12 .* n_d12' - n_d .* n_d'
            inter2 = n_d12 .* n_d12' + n_d .* n_d'
            # Intrapair Electron Correlation
            for l in 1:pp.ndoc
                ldx = pp.no1 + l
                # inicio y fin de los orbitales acoplados a los fuertemente ocupados
                ll = pp.no1 + pp.ndns + pp.ncwo*(pp.ndoc - l) + 1
                ul = pp.no1 + pp.ndns + pp.ncwo*(pp.ndoc - l + 1)

                inter[ldx,ldx] = 0
                inter[ldx,ll:ul] .= 0
                inter[ll:ul,ldx] .= 0
                inter[ll:ul,ll:ul] .= 0
                inter2[ldx,ldx] = 0
                inter2[ldx,ll:ul] .= 0
                inter2[ll:ul,ldx] .= 0
                inter2[ll:ul,ll:ul] .= 0
            end

	    inter[pp.nbeta+1:end,pp.nbeta+1:end] .= 0
            inter[1:pp.nalpha,1:pp.nalpha] .= 0
	    inter2[pp.nbeta+1:end,pp.nbeta+1:end] .= 0
            inter2[1:pp.nalpha,:] .= 0
            inter2[:,1:pp.nalpha] .= 0

            Pi_d = inter - inter2

	    @tullio Dab[p,q,r,t] += -Pi_d[p,r]*Id[p,q]*Id[r,t]
        end
    end

    return Daa, Dab

end
