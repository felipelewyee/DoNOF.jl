module pnof

using Tullio
using TensorOperations
using LinearAlgebra

function PNOFi_selector(n,p)
    #if(p.ipnof==5):
    #    cj12,ck12 = CJCKD5(n,p)
    if p.ipnof==7
        cj12,ck12 = CJCKD7(n,p.ista,p.no1,p.ndoc,p.nsoc,p.nbeta,p.nalpha,p.ndns,p.ncwo,p.MSpin)
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
    #cj12 = 2*np.outer(n,n)
    #ck12 = np.outer(n,n) + np.outer(fi,fi)

    # Intrapair Electron Correlation

    if MSpin==0 && nsoc>1
	n_beta_alpha = view(n,nbeta+1:nalpha)
	ck12_beta_alpha = view(ck,nbeta+1:nalpha,nbeta+1:nalpha)
	@tullio ck12_beta_alpha[i,j] := 2*n_beta_alpha[i]*n_beta_alpha[j]
        #ck12[nbeta:nalpha,nbeta:nalpha] = 2*np.outer(n[nbeta:nalpha],n[nbeta:nalpha])
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
        #ck12[ll:ul,ll:ul] = -sqrt.(outer(n[ll:ul],n[ll:ul]))
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
                dfi_dgamma[i,k] = 1/(2*a)*(1-2*n[i])*dn_dgamma[i][k]
            else
                dfi_dgamma[i,k] = 2*(1-2*n[i])*dn_dgamma[i][k]
            end
	end
    end

    # Interpair Electron correlation #

    #Dcj12r = np.zeros((nbf5,nbf5,nv))
    #Dck12r = np.zeros((nbf5,nbf5,nv))
    @tullio Dj12r[i,j,k] := 2*dn_dgamma[i,k]*n[j]
    @tullio Dk12r[i,j,k] := dn_dgamma[i,k]*n[j] + dfi_dgamma[i,k]*fi[j]
    #for k in 1:nv
    #    Dcj12r[:,:,k] = 2*np.outer(dn_dgamma[:,k],n)
    #    Dck12r[:,:,k] = np.outer(dn_dgamma[:,k],n) + np.outer(dfi_dgamma[:,k],fi)
    #Dcj12r = 2*np.einsum('ik,j->ijk',dn_dgamma,n)    
    #Dck12r = np.einsum('ik,j->ijk',dn_dgamma,n) + np.einsum('ik,j->ijk',dfi_dgamma,fi)    

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

        for k in 1:nv
            Dck12r[ldx,ll:ul,k] = 1/2 * 1/np.sqrt(a) * dn_dgamma[ldx,k] * np.sqrt(n[ll:ul])
            Dck12r[ll:ul,ldx,k] = 1/2 * 1/np.sqrt(b) * dn_dgamma[ll:ul,k] * np.sqrt(n[ldx])
            Dck12r[ll:ul,ll:ul,k] = - 1/2 * np.outer(1/np.sqrt(b)*dn_dgamma[ll:ul,k],np.sqrt(n[ll:ul]))
        end	
        #Dck12r[ldx,ll:ul,:nv] = 1/2 * 1/np.sqrt(a) * np.einsum('j,i->ij',dn_dgamma[ldx,:nv],np.sqrt(n[ll:ul]))
        #Dck12r[ll:ul,ldx,:nv] = 1/2 * np.einsum('i,ij->ij', 1/np.sqrt(b),dn_dgamma[ll:ul,:nv])*np.sqrt(n[ldx])

        #for k in range(nv):
        #    Dck12r[ll:ul,ll:ul,k] = - 1/2 * np.einsum('i,i,j->ij',1/np.sqrt(b),dn_dgamma[ll:ul,k],np.sqrt(n[ll:ul]))
    end
    return Dcj12r,Dck12r

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
                n[ll+iw] .*= sin(gamma[ndoc+(ncwo-1)*(i-1)+iw])^2
	        n[ll+iw:ul] .*= cos(gamma[ndoc+(ncwo-1)*(i-1)+iw])^2
	    end
	end

        for i in 1:ndoc
            # dn_g/dgamma_g
            dn_dgamma[no1+i,i] = dni_dgammai[no1+i]

            # dn_pi/dgamma_g
	    ll = no1 + ndns + ncwo*(ndoc - i) + 1
            ul = no1 + ndns + ncwo*(ndoc - i + 1)
            dn_dgamma[ll:ul,i] = -dni_dgammai[no1+i]
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
			    dn_dgamma[ip][ndoc+(ncwo-1)*(i-1)+iw] *= np.cos(gamma[ndoc+(ncwo-1)*(i-1)+jw])^2
                        end
                    if ip-ll<ncwo-1
                        dn_dgamma[ip][ndoc+(ncwo-1)*(i+1)+iw] *= np.sin(gamma[ndoc+(ncwo-1)*(i+1)+(ip-ll)])^2
                    end
                    end
                end
            end

            # dn_pi/dgamma_i
            for iw in 1:ncwo-1
                dn_dgamma[ll+iw-1][ndoc+(ncwo-1)*(i-1)+iw] = 1 - n[no1+i]
                for jw in 1:iw+1
                    if jw==iw
                        dn_dgamma[ll+iw][ndoc+(ncwo-1)*(i-1)+iw] *= np.sin(2*gamma[ndoc+(ncwo-1)*(i-1)+jw])
                    else
                        dn_dgamma[ll+iw][ndoc+(ncwo-1)*(i-1)+iw] *= np.cos(gamma[ndoc+(ncwo-1)*(i-1)+jw])^2
		    end
		end
	    end
        end


        end
    return n,dn_dgamma


end

function calce(gamma,J_MO,K_MO,H_core,p)

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin)
    cj12,ck12 = PNOFi_selector(n,p)

    E = 0

    if p.MSpin==0

        # 2H + J
	n_beta = view(n,1:p.nbeta)
	n_alpha = view(n,p.nbeta+1:p.nalpha)
	n_nbf5 = view(n,p.nalpha+1:p.nbf5)
	H_beta = view(H_core,1:p.nbeta)
	H_alpha= view(H_core,p.nbeta+1:p.nalpha)
	H_nbf5 = view(H_core,p.nalpha+1:p.nbf5)
	J_MO_beta = view(J_MO,1:p.nbeta,1:p.nbeta)
	J_MO_nbf5 = view(J_MO,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
	@tullio E += n_beta[i]*(2*H_beta[i] + J_MO_beta[i,i])
	#E = E + np.einsum('i,i',n[:p.nbeta],2*H_core[:p.nbeta]+np.diagonal(J_MO)[:p.nbeta],optimize=True) # [0,Nbeta]
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

#    elif(not p.MSpin==0):
#        E = 0
#
#        # 2H + J
#        E = E + np.einsum('i,i',n[:p.nbeta],2*H_core[:p.nbeta]+np.diagonal(J_MO)[:p.nbeta],optimize=True) # [0,Nbeta]
#        E = E + np.einsum('i,i',n[p.nbeta:p.nalpha],2*H_core[p.nbeta:p.nalpha],optimize=True)               # (Nbeta,Nalpha]
#        E = E + np.einsum('i,i',n[p.nalpha:p.nbf5],2*H_core[p.nalpha:p.nbf5]+np.diagonal(J_MO)[p.nalpha:p.nbf5],optimize=True) # (Nalpha,Nbf5)
#
#        #C^J JMO
#        np.fill_diagonal(cj12,0) # Remove diag.
#        E = E + np.einsum('ij,ji->',cj12[:p.nbeta,:p.nbeta],J_MO[:p.nbeta,:p.nbeta],optimize=True) # sum_ij
#        E = E + np.einsum('ij,ji->',cj12[:p.nbeta,p.nalpha:p.nbf5],J_MO[p.nalpha:p.nbf5,:p.nbeta],optimize=True) # sum_ij
#        E = E + np.einsum('ij,ji->',cj12[p.nalpha:p.nbf5,:p.nbeta],J_MO[:p.nbeta,p.nalpha:p.nbf5],optimize=True) # sum_ij
#        E = E + np.einsum('ij,ji->',cj12[p.nalpha:p.nbf5,p.nalpha:p.nbf5],J_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5],optimize=True) # sum_ij
#
#        #C^K KMO
#        np.fill_diagonal(ck12,0) # Remove diag.
#        E = E - np.einsum('ij,ji->',ck12[:p.nbeta,:p.nbeta],K_MO[:p.nbeta,:p.nbeta],optimize=True) # sum_ij
#        E = E - np.einsum('ij,ji->',ck12[:p.nbeta,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,:p.nbeta],optimize=True) # sum_ij
#        E = E - np.einsum('ij,ji->',ck12[p.nalpha:p.nbf5,:p.nbeta],K_MO[:p.nbeta,p.nalpha:p.nbf5],optimize=True) # sum_ij
#        E = E - np.einsum('ij,ji->',ck12[p.nalpha:p.nbf5,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5],optimize=True) # sum_ij
#
#        #n JMO
#        E = E + 2*np.einsum('i,ji->',n[:p.nbeta],J_MO[p.nbeta:p.nalpha,:p.nbeta],optimize=True) # sum_ij
#        E = E + 2*np.einsum('i,ji->',n[p.nalpha:p.nbf5],J_MO[p.nbeta:p.nalpha,p.nalpha:p.nbf5],optimize=True) # sum_ij
#        E = E + 0.5*(np.einsum('i,ji->',n[p.nbeta:p.nalpha],J_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha],optimize=True) - np.einsum('i,ii->',n[p.nbeta:p.nalpha],J_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha],optimize=True))
#
#        #n KMO
#        E = E - np.einsum('i,ji->',n[:p.nbeta],K_MO[p.nbeta:p.nalpha,:p.nbeta],optimize=True) # sum_ij
#        E = E - np.einsum('i,ji->',n[p.nalpha:p.nbf5],K_MO[p.nbeta:p.nalpha,p.nalpha:p.nbf5],optimize=True) # sum_ij
#        E = E - np.einsum('i,ji->',n[p.nbeta:p.nalpha],K_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha],optimize=True) - np.einsum('i,ii->',n[p.nbeta:p.nalpha],K_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha],optimize=True)
    end
    return E
end

end
