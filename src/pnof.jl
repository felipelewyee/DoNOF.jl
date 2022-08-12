function CJCKD_muller(n)
    @tullio cj12[i,j] := 2*n[i]*n[j]
    @tullio ck12[i,j] := sqrt(n[i]*n[j]) #revisar sintaxis!!!
    return cj12,ck12
end

function der_CJCKD_muller(n,dn_dgamma)
    Dcj12r[i,j,k] = dn_dgamma[i,k]*n[j]
    Dck12r[i,j,k] = dn_dgamma[i,k]*sqrt.(n[j])

    return Dcj12r,Dck12r
end

function PNOFi_selector(n,p)
    #if(p.ipnof==5):
    #    cj12,ck12 = CJCKD5(n,p)
    if p.ipnof==7
        cj12,ck12 = CJCKD7(n,p.ista,p.no1,p.ndoc,p.nsoc,p.nbeta,p.nalpha,p.ndns,p.ncwo,p.MSpin)
    elseif p.ipnof==8
        cj12,ck12 = CJCKD8(n,p.ista,p.no1,p.ndoc,p.nsoc,p.nbeta,p.nalpha,p.ndns,p.ncwo,p.MSpin,p.lamb)
    end

    return cj12,ck12
end

function der_PNOFi_selector(n,dn_dgamma,p)
    #if p.ipnof==5
    #    Dcj12r,Dck12r = der_CJCKD5(n,dn_dgamma,p)
    if p.ipnof==7
        Dcj12r,Dck12r = der_CJCKD7(n,p.ista,dn_dgamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo)
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

function CJCKD8(n,ista,no1,ndoc,nsoc,nbeta,nalpha,ndns,ncwo,MSpin,lamb)

    nbf5 = size(n)[1]
    delta = zeros(nbf5,nbf5)
    for i in 1:nbf5
        for j in 1:nbf5
            delta[i,j] = lamb*min(n[i]*n[j],(1-n[i])*(1-n[j]))
        end
    end
    ppi = zeros(nbf5,nbf5)
    for i in 1:nbf5
        for j in 1:nbf5
            ppi[i,j] = sqrt((n[i]*(1-n[j])+delta[i,j]) * (n[j]*(1-n[i])+delta[j,i]))
        end
    end

    # Interpair Electron correlation #

    @tullio cj12[i,j] := 2*n[i]*n[j] - 2*delta[i,j]
    @tullio ck12[i,j] := n[i]*n[j] - delta[i,j] + ppi[i,j]

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

function der_CJCKD8(n,ista,dn_dgamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo)

    nbf5 = size(n)[1]
    delta = zeros(nbf5,nbf5)
    for i in 1:nbf5
        for j in 1:nbf5
            delta[i,j] = lamb*min(n[i]*n[j],(1-n[i])*(1-n[j]))
        end
    end
    ppi = zeros(nbf5,nbf5)
    for i in 1:nbf5
        for j in 1:nbf5
            ppi[i,j] = sqrt((n[i]*(1-n[j])+delta[i,j]) * (n[j]*(1-n[i])+delta[j,i]))
        end
    end

    # Interpair Electron correlation #

    @tullio cj12[i,j] := 2*n[i]*n[j] - delta[i,j]
    @tullio ck12[i,j] := n[i]*n[j] - delta[i,j] + ppi[i,j]

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

function f(mu,gamma,nalpha,no1)

    N=0.0
    for gi in gamma
        N += (erf(gi+mu[1]) + 1)/2
    end

    return (nalpha-N-no1)^2

end

function compute_mu(gamma,nalpha,no1)
    mu = [0.0]
    res = optimize(mu->f(mu,gamma,nalpha,no1),mu)
    #mu = res.x
    mu = Optim.minimizer(res)

   return mu[1]

end


function ocupacion_muller(nv,nbf,no1,nalpha,gamma,p)
    n = zeros(nbf)
    n[1:no1] .= 1
    if(p.EBI)
        mu = compute_mu(gamma,nalpha,no1)

        n[no1+1:no1+nv] = (erf.(gamma .+ mu) .+ 1)./2
        dn_dgamma = zeros(nbf,nv)
        dn_dgamma = zeros(nbf,nv)
    else
        dni_dgammai = zeros(nbf)

        n[no1+1:no1+nv] = cos.(gamma[1:nv]).^2
        dni_dgammai = - sin.(2*gamma)

        dn_dgamma = zeros.(nbf,nv)

        dn_dgamma = Diagonal(dni_dgammai)
    end
    return n,dn_dgamma

end


function ocupacion(gamma,no1,ndoc,nalpha,nv,nbf5,ndns,ncwo,HighSpin)

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

    if(ncwo==1)
        dn_dgamma = zeros(nbf5,nv)

        for i in 1:ndoc
            dn_dgamma[no1+i,i] = dni_dgammai[no1+i]
            #cwo
            icf = nalpha + ndoc - i + 1
	    n[icf] = 1/2*(sin(gamma[i]))^2
            dni_dgammai[icf]  = 1/2*sin(2*gamma[i])
            dn_dgamma[icf,i] = dni_dgammai[icf]
        end
    else
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


        end
    return n,dn_dgamma


end


function calce_muller(gamma,J_MO,K_MO,H_core,nv,nbf,no1,nalpha,p)
    n,dn_dgamma = ocupacion_muller(nv,nbf,no1,nalpha,gamma,p)
    cj12,ck12 = CJCKD_muller(n)
    E = 0
    #nH
    @tullio E += n[i]*2*H_core[i]
    #C^J JMO
    @tullio E += cj12[i,j]*J_MO[j,l]
    #C^K KMO
    @tullio E += ck12[i,j]*K_MO[j,l]

    return E
end


function calce(gamma,J_MO,K_MO,H_core,p)

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin)
    cj12,ck12 = PNOFi_selector(n,p)

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
	#E = E + np.einsum('i,i',n[:p.nbeta],2*H_core[:p.nbeta]+np.diagonal(J_MO)[:p.nbeta],optimize=True) # [0,Nbeta]
        @tullio E += n_beta[i]*(2*H_beta[i] + J_MO_beta[i,i])
	@tullio E += n_alpha[i]*2*H_alpha[i]
        #E = E + np.einsum('i,i',n[p.nbeta:p.nalpha],2*H_core[p.nbeta:p.nalpha],optimize=True)               # (Nbeta,Nalpha]
	@tullio E += n_nbf5[i]*(2*H_nbf5[i]+J_MO_nbf5[i,i])
        #E = E + np.einsum('i,i',n[p.nalpha:p.nbf5],2*H_core[p.nalpha:p.nbf5]+np.diagonal(J_MO)[p.nalpha:p.nbf5],optimize=True) # (Nalpha,Nbf5)

        #C^J JMO
        cj12[diagind(cj12)] .= 0
	@tullio E += cj12[i,j]*J_MO[j,i]
        #np.fill_diagonal(cj12,0) # Remove diag.
        #E = E + np.einsum('ij,ji->',cj12,J_MO,optimize=True) # sum_ij

        #C^K KMO
        ck12[diagind(ck12)] .= 0            
	@tullio E += -ck12[i,j]*K_MO[j,i]
        #np.fill_diagonal(ck12,0) # Remove diag.
        #E = E - np.einsum('ij,ji->',ck12,K_MO,optimize=True) # sum_ij

    elseif p.MSpin!=0

        # 2H + J
        @tullio E += n_beta[i]*(2*H_beta[i] + J_MO_beta[i,i])
#        E = E + np.einsum('i,i',n[:p.nbeta],2*H_core[:p.nbeta]+np.diagonal(J_MO)[:p.nbeta],optimize=True) # [0,Nbeta]
	@tullio E += n_alpha[i]*2*H_alpha[i]
#        E = E + np.einsum('i,i',n[p.nbeta:p.nalpha],2*H_core[p.nbeta:p.nalpha],optimize=True)               # (Nbeta,Nalpha]
	@tullio E += n_nbf5[i]*(2*H_nbf5[i]+J_MO_nbf5[i,i])
#        E = E + np.einsum('i,i',n[p.nalpha:p.nbf5],2*H_core[p.nalpha:p.nbf5]+np.diagonal(J_MO)[p.nalpha:p.nbf5],optimize=True) # (Nalpha,Nbf5)
#
#        #C^J JMO
#        np.fill_diagonal(cj12,0) # Remove diag.
        cj12[diagind(cj12)] .= 0
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
#        E = E + np.einsum('ij,ji->',cj12[:p.nbeta,:p.nbeta],J_MO[:p.nbeta,:p.nbeta],optimize=True) # sum_ij
	@tullio E += cj12_beta_nbf5[i,j]*J_MO_nbf5_beta[j,i]
#        E = E + np.einsum('ij,ji->',cj12[:p.nbeta,p.nalpha:p.nbf5],J_MO[p.nalpha:p.nbf5,:p.nbeta],optimize=True) # sum_ij
	@tullio E += cj12_nbf5_beta[i,j]*J_MO_beta_nbf5[j,i]
#        E = E + np.einsum('ij,ji->',cj12[p.nalpha:p.nbf5,:p.nbeta],J_MO[:p.nbeta,p.nalpha:p.nbf5],optimize=True) # sum_ij
	@tullio E += cj12_nbf5_nbf5[i,j]*J_MO_nbf5_nbf5[j,i]
#        E = E + np.einsum('ij,ji->',cj12[p.nalpha:p.nbf5,p.nalpha:p.nbf5],J_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5],optimize=True) # sum_ij
#
#        #C^K KMO
#        np.fill_diagonal(ck12,0) # Remove diag.
        ck12[diagind(ck12)] .= 0
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
#        E = E - np.einsum('ij,ji->',ck12[:p.nbeta,:p.nbeta],K_MO[:p.nbeta,:p.nbeta],optimize=True) # sum_ij
	@tullio E += -ck12_beta_nbf5[i,j]*K_MO_nbf5_beta[j,i]
#        E = E - np.einsum('ij,ji->',ck12[:p.nbeta,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,:p.nbeta],optimize=True) # sum_ij
	@tullio E += -ck12_nbf5_beta[i,j]*K_MO_beta_nbf5[j,i]
#        E = E - np.einsum('ij,ji->',ck12[p.nalpha:p.nbf5,:p.nbeta],K_MO[:p.nbeta,p.nalpha:p.nbf5],optimize=True) # sum_ij
	@tullio E += -ck12_nbf5_nbf5[i,j]*K_MO_nbf5_nbf5[j,i]
#        E = E - np.einsum('ij,ji->',ck12[p.nalpha:p.nbf5,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5],optimize=True) # sum_ij
#
#        #n JMO
	@tullio E += 2*n_beta[i]*J_MO_alpha_beta[j,i]
#        E = E + 2*np.einsum('i,ji->',n[:p.nbeta],J_MO[p.nbeta:p.nalpha,:p.nbeta],optimize=True) # sum_ij
	@tullio E += 2*n_nbf5[i]*J_MO_alpha_nbf5[j,i]
#        E = E + 2*np.einsum('i,ji->',n[p.nalpha:p.nbf5],J_MO[p.nbeta:p.nalpha,p.nalpha:p.nbf5],optimize=True) # sum_ij
	@tullio E += 0.5*n_alpha[i]*J_MO_alpha_alpha[j,i]
	@tullio E += -0.5*n_alpha[i]*J_MO_alpha_alpha[i,i]
#        E = E + 0.5*(np.einsum('i,ji->',n[p.nbeta:p.nalpha],J_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha],optimize=True) - np.einsum('i,ii->',n[p.nbeta:p.nalpha],J_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha],optimize=True))
#
#        #n KMO
	@tullio E += -n_beta[i]*K_MO_alpha_beta[j,i]
#        E = E - np.einsum('i,ji->',n[:p.nbeta],K_MO[p.nbeta:p.nalpha,:p.nbeta],optimize=True) # sum_ij
	@tullio E += -n_nbf5[i]*K_MO_alpha_nbf5[j,i]
#        E = E - np.einsum('i,ji->',n[p.nalpha:p.nbf5],K_MO[p.nbeta:p.nalpha,p.nalpha:p.nbf5],optimize=True) # sum_ij
	@tullio E += -n_alpha[i]*K_MO_alpha_alpha[j,i]
	@tullio E += n_alpha[i]*K_MO_alpha_alpha[i,i]
#        E = E - np.einsum('i,ji->',n[p.nbeta:p.nalpha],K_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha],optimize=True) - np.einsum('i,ii->',n[p.nbeta:p.nalpha],K_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha],optimize=True)
    end
    #for i in 1:size(dn_dgamma)[1]
    #    for j in 1:size(dn_dgamma)[2]
    #        println(i," ",j," ",dn_dgamma[i,j])
    #    end
    #end 
    #stop
    return E
end


function calcg_muller(gamma,J_MO,K_MO,H_core,nv,nbf,no1)
    grad = zeros(nv)
    n,dn_dgamma = ocupacion_muller(nv,nbf,no1,nalpha,gamma,p)
    Dcj12r,Dck12r = der_CJCKD_muller(n,dn_dgamma)


    #nH
    grad[k] += dn_dgamma[i,k]*2*H_core[i]
    #C^J JMO
    grad[k] += Dcj12r[i,j,k]*J_MO[j,i]
    #C^K KMO
    grad[k] += Dck12r[i,j,k]*K_MO[j,i]

   return grad
end


function calcg(gamma,J_MO,K_MO,H_core,p)

    grad = zeros(p.nv)

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin)
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
        #diag = np.diag_indices(p.nbf5)
        #Dcj12r[diag] = 0
	for i in 1:p.nv
	    for j in 1:p.nbf5
	        Dcj12r[j,j,i] = 0
            end 
        end
	Dcj12r_beta = view(Dcj12r,p.no1+1:p.nbeta,1:p.nbf5,1:p.nv)
	Dcj12r_nbf5 = view(Dcj12r,p.nalpha+1:p.nbf5,1:p.nbf5,1:p.nv)
	J_MO_beta = view(J_MO,1:p.nbf5,p.no1+1:p.nbeta)
	J_MO_nbf5 = view(J_MO,1:p.nbf5,p.nalpha+1:p.nbf5)
	@tullio grad[k] += 2*Dcj12r_beta[i,j,k]*J_MO_beta[j,i]
        #grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.no1:p.nbeta,:p.nbf5,:p.nv],J_MO[:p.nbf5,p.no1:p.nbeta],optimize=True)

	@tullio grad[k] += 2*Dcj12r_nbf5[i,j,k]*J_MO_nbf5[j,i]
        #grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.nalpha:p.nbf5,:p.nbf5,:p.nv],J_MO[:p.nbf5,p.nalpha:p.nbf5],optimize=True)

        # -2 dCK_dgamma K_MO
        #diag = np.diag_indices(p.nbf5)
	for i in 1:p.nv
	    for j in 1:p.nbf5
	        Dck12r[j,j,i] = 0
            end
        end
	Dck12r_beta = view(Dck12r,p.no1+1:p.nbeta,1:p.nbf5,1:p.nv)
	Dck12r_nbf5 = view(Dck12r,p.nalpha+1:p.nbf5,1:p.nbf5,1:p.nv)
	K_MO_beta = view(K_MO,1:p.nbf5,p.no1+1:p.nbeta)
	K_MO_nbf5 = view(K_MO,1:p.nbf5,p.nalpha+1:p.nbf5)
        #Dck12r[diag] = 0
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
        for i in 1:p.nv
            for j in 1:p.nbf5
                Dcj12r[j,j,i] = 0
            end
        end
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
#        Dcj12r[np.diag_indices(p.nbf5)] = 0
#        grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.no1:p.nbeta,:p.nbeta,:p.nv],J_MO[:p.nbeta,p.no1:p.nbeta],optimize=True)
#        grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.no1:p.nbeta,p.nalpha:p.nbf5,:p.nv],J_MO[p.nalpha:p.nbf5,p.no1:p.nbeta],optimize=True)
#        grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.nalpha:p.nbf5,:p.nbeta,:p.nv],J_MO[:p.nbeta,p.nalpha:p.nbf5],optimize=True)
#        grad += 2*np.einsum('ijk,ji->k',Dcj12r[p.nalpha:p.nbf5,p.nalpha:p.nbf5,:p.nv],J_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5],optimize=True)
#
#        # -2 dCK_dgamma K_MO
        for i in 1:p.nv
            for j in 1:p.nbf5
                Dck12r[j,j,i] = 0
            end
        end
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
#        Dck12r[np.diag_indices(p.nbf5)] = 0
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

function calcorbe(y,gamma,C,H,I,b_mnl,p)

    Cnew = rotate_orbital(y,C,p)
    print("Cnew calculado")
    J_MO,K_MO,H_core = computeJKH_MO(Cnew,H,I,b_mnl,p)
    print("integrales MO calculadas")
    if(p.muller)
	E = calce_muller(gamma,J_MO,K_MO,H_core,nv,nbf,no1,nalpha,p)
    else
        E = calce(gamma,J_MO,K_MO,H_core,p)
    print("Energía calculada")
    end
    return E

end

function calcorbg(y,gamma,C,H,I,b_mnl,pa)

    Cnew = rotate_orbital(y,C,pa)
    print("Cnew calculado")
    if(pa.gpu)
        Cnew = CuArray(Cnew)
	print("CNew con CUDA definido")
        H = CuArray(H)
	print("H con CUDA definido")
    end
    if(pa.muller)
	n,dn_dgamma = ocupacion_muller(pa.nv,pa.nbf,pa.no1,pa.nalpha,gamma,pa) 
	print("n, dngamma calculados")
	cj12,ck12 = CJCKD_muller(n)
	print("matrices cj12 y ck12 calculadas")
    else
        n,dn_dgamma = ocupacion(gamma,pa.no1,pa.ndoc,pa.nalpha,pa.nv,pa.nbf5,pa.ndns,pa.ncwo,pa.HighSpin)
        cj12,ck12 = PNOFi_selector(n,pa)
    end

    @tullio H_mat[i,j] := Cnew[m,i] * H[m,nn] * Cnew[nn,j]
    print("H_mar calculada")
    if pa.RI
        @tullio b_MO[p,q,l] := Cnew[m,p]*Cnew[nn,q]*b_mnl[m,nn,l]
        print("si se imprime esto algo salió mal") 
    else
    	@tullio I_MO[p,q,r,t] := Cnew[m,p]*Cnew[nn,q]*I[m,nn,s,l]*Cnew[s,r]*Cnew[l,t]
	print("I_MO calculada")
    end
   
    grad = zeros(pa.nbf,pa.nbf)
    print("se definió grad") 
    #if(!pa.muller)
    #	print("entró a pa.muller")
    #    cj12[diagind(cj12)] .= 0
    #    ck12[diagind(ck12)] .= 0
    #	print("si sale esto es porque algo salió mal")
    #end 

    if pa.gpu
        grad = CuArray(grad)
        cj12 = CuArray(cj12)
        ck12 = CuArray(ck12)
        n = CuArray(n)
    print("grad,cj12,ck12 calculados para gpu")
    end 

    n_beta =        view(n,1:pa.nbeta)
    n_alpha =       view(n,pa.nalpha+1:pa.nbf5)
    Hmat_nbf5 =      view(H_mat,1:pa.nbf,1:pa.nbf5)
    Hmat_nbf5t =     view(H_mat,1:pa.nbf5,1:pa.nbf)
    grad_nbf5 =      view(grad,1:pa.nbf,1:pa.nbf5)
    grad_nbf5t =     view(grad,1:pa.nbf5,1:pa.nbf)
    grad_nbeta =     view(grad,1:pa.nbf,1:pa.nbeta)
    grad_nbetat =    view(grad,1:pa.nbeta,1:pa.nbf)
    grad_nalpha =    view(grad,1:pa.nbf,pa.nalpha+1:pa.nbf5)
    grad_nalphat =   view(grad,pa.nalpha+1:pa.nbf5,1:pa.nbf)
    if pa.RI
    b_nbf_beta =     view(b_MO,1:pa.nbf,1:pa.nbeta,1:pa.nbfaux)
    b_nbf_alpha =    view(b_MO,1:pa.nbf,pa.nalpha+1:pa.nbf5,1:pa.nbfaux)
    b_nbeta_beta =   view(b_MO,1:pa.nbeta,1:pa.nbeta,1:pa.nbfaux)
    b_nalpha_alpha = view(b_MO,pa.nalpha+1:pa.nbf5,pa.nalpha+1:pa.nbf5,1:pa.nbfaux)
    b_nbf_nbf5 =     view(b_MO,1:pa.nbf,1:pa.nbf5,1:pa.nbfaux)
    b_nbf5_nbf5 =    view(b_MO,1:pa.nbf5,1:pa.nbf5,1:pa.nbfaux)
    else
    I_nb_nb_nb =    view(I_MO,1:pa.nbf,1:pa.nbeta,1:pa.nbeta,1:pa.nbeta)
    I_na_na_na =    view(I_MO,1:pa.nbf,pa.nalpha+1:pa.nbf5,pa.nalpha+1:pa.nbf5,pa.nalpha+1:pa.nbf5)
    I_nbf5_nbf5_nbf5 =    view(I_MO,1:pa.nbf,1:pa.nbf5,1:pa.nbf5,1:pa.nbf5)
    end
    print("definiciones para lo que no es Muller calculadas")
    if pa.RI
        if(pa.MSpin==0)
	    if(pa.muller)
		# 2ndH/dy_ab
	        @tullio grad[a,b] += +4*n[b]*H_mat[a,b]
                @tullio grad[a,b] += -4*n[a]*H_mat[a,b]
		# C^J_pq dJ_pq/dy_ab 
                @tullio grad += +4*cj12[b,q]*b_MO[a,b,k]*b_MO[q,q,k]
                @tullio grad += -4*cj12[a,q]*b_MO[b,a,k]*b_MO[q,q,k]

                # -C^K_pq dK_pq/dy_ab 
                @tullio grad += -4*ck12[b,q]*b_MO[a,q,k]*b_MO[b,q,k]
                @tullio grad += +4*ck12[a,q]*b_MO[b,q,k]*b_MO[a,q,k]

	    else
                # 2ndH/dy_ab
                @tullio grad_nbf5[a,b]  += +4*n[b]*Hmat_nbf5[a,b] 
                @tullio grad_nbf5t[a,b] += -4*n[a]*Hmat_nbf5t[a,b] 

                # dJ_pp/dy_ab
                @tullio grad_nbeta[a,b] += +4*n_beta[b]*b_nbf_beta[a,b,k]*b_nbeta_beta[b,b,k]
                @tullio grad_nalpha[a,b] += +4*n_alpha[b]*b_nbf_alpha[a,b,k]*b_nalpha_alpha[b,b,k]
	        @tullio grad_nbetat[a,b] += -4*n_beta[a]*b_nbf_beta[b,a,k]*b_nbeta_beta[a,a,k]
                @tullio grad_nalphat[a,b] += -4*n_alpha[a]*b_nbf_alpha[b,a,k]*b_nalpha_alpha[a,a,k]

                # C^J_pq dJ_pq/dy_ab 
                @tullio grad_nbf5[a,b] += +4*cj12[b,q]*b_nbf_nbf5[a,b,k]*b_nbf5_nbf5[q,q,k]
                @tullio grad_nbf5t[a,b] += -4*cj12[a,q]*b_nbf_nbf5[b,a,k]*b_nbf5_nbf5[q,q,k]

                # -C^K_pq dK_pq/dy_ab 
                @tullio grad_nbf5[a,b] += -4*ck12[b,q]*b_nbf_nbf5[a,q,k]*b_nbf5_nbf5[b,q,k]
                @tullio grad_nbf5t[a,b] += +4*ck12[a,q]*b_nbf_nbf5[b,q,k]*b_nbf5_nbf5[a,q,k]
            end
	end
    else
        if(pa.MSpin==0)
            if(pa.muller)
	        # 2ndH/dy_ab
                @tullio grad[a,b]  += +4*n[b]*Hmat_nbf5[a,b]
                @tullio grad[a,b] += -4*n[a]*Hmat_nbf5t[a,b]
		print("2ndH/dy_ab calculadas")

                # C^J_pq dJ_pq/dy_ab 
                @tullio grad[a,b] += +4*cj12[b,q]*I_MO[a,b,q,q]
                @tullio grad[a,b] += -4*cj12[a,q]*I_MO[b,a,q,q]
                print("Relacionadas con J calculadas")
                # -C^K_pq dK_pq/dy_ab 
                @tullio grad[a,b] += -4*ck12[b,q]*I_MO[a,q,b,q]
                @tullio grad[a,b] += +4*ck12[a,q]*I_MO[b,q,a,q]
                print("Relacionadas con K calculadas")
	    else
                # 2ndH/dy_ab
                @tullio grad_nbf5[a,b]  += +4*n[b]*Hmat_nbf5[a,b] 
                @tullio grad_nbf5t[a,b] += -4*n[a]*Hmat_nbf5t[a,b] 

                # dJ_pp/dy_ab
	        @tullio grad_nbeta[a,b] += +4*n_beta[b]*I_nb_nb_nb[a,b,b,b]
                @tullio grad_nalpha[a,b] += +4*n_alpha[b]*I_na_na_na[a,b,b,b]
                @tullio grad_nbetat[a,b] += -4*n_beta[a]*I_nb_nb_nb[b,a,a,a]
                @tullio grad_nalphat[a,b] += -4*n_alpha[a]*I_na_na_na[b,a,a,a]

                # C^J_pq dJ_pq/dy_ab
	        @tullio grad_nbf5[a,b] += +4*cj12[b,q]*I_nbf5_nbf5_nbf5[a,b,q,q]
	        @tullio grad_nbf5t[a,b] += -4*cj12[a,q]*I_nbf5_nbf5_nbf5[b,a,q,q]

                # -C^K_pq dK_pq/dy_ab 
       	        @tullio grad_nbf5[a,b] += -4*ck12[b,q]*I_nbf5_nbf5_nbf5[a,q,b,q]
                @tullio grad_nbf5t[a,b] += +4*ck12[a,q]*I_nbf5_nbf5_nbf5[b,q,a,q]
		print("si aparece esto es porque algo va mal")
            end
        end
    end

    grad = Array(grad)
    grads = zeros(pa.nvar)
    n = 1

    for i in 1:pa.nbf5
        for j in i+1:pa.nbf
            grads[n] = grad[i,j]
            n += 1
         end
    end
    
    grad = grads
    print("grad calculada")
    return grad

end
