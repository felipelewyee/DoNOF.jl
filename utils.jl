module utils

using Tullio
using LinearAlgebra

include("integrals.jl")

function ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)

    #J,K = integrals.computeJKj(C,I,b_mnl,p)
    J,K = integrals.JKj_Full(C,I,b_mnl,p)
    
    if p.MSpin==0
        #F = computeF_RC_driver(J,K,n,H,cj12,ck12,p)
        F = computeF_RC_CPU(J,K,n,H,cj12,ck12,p)
#    elseif not p.MSpin==0
#        F = computeF_RO_driver(J,K,n,H,cj12,ck12,p)
    end

    elag = computeLagrange(F,C,p)

    E = computeE_elec(H,C,n,elag,p)

    sumdiff,maxdiff = computeLagrangeConvergency(elag)

    return E,elag,sumdiff,maxdiff


end

function computeF_RC_CPU(J,K,n,H,cj12,ck12,p)

    # Matriz de Fock Generalizada
    F = zeros(p.nbf5,p.nbf,p.nbf)

    ini = 1
    if(p.no1>1)
        ini = p.no1 + 1
    end

    # nH
    @tullio F[i,m,s] += n[i]*H[m,s]
    #F += np.einsum('i,mn->imn',n,H,optimize=True)        # i = [1,nbf5]

    # nJ
    F_ini_beta = view(F,ini:p.nbeta,:,:)
    n_ini_beta = view(n,ini:p.nbeta)
    J_ini_beta = view(J,ini:p.nbeta,:,:)
    @tullio F_ini_beta[i,m,n] += n_ini_beta[i]*J_ini_beta[i,m,n]
    #F[ini:p.nbeta,:,:] += np.einsum('i,imn->imn',n[ini:p.nbeta],J[ini:p.nbeta,:,:],optimize=True)        # i = [ini,nbeta]
    F_alpha_nbf5 = view(F,p.nalpha+1:p.nbf5,:,:)
    n_alpha_nbf5 = view(n,p.nalpha+1:p.nbf5)
    J_alpha_nbf5 = view(J,p.nalpha+1:p.nbf5,:,:)
    @tullio F_alpha_nbf5[i,m,n] += n_alpha_nbf5[i]*J_alpha_nbf5[i,m,n]
    #F[p.nalpha:p.nbf5,:,:] += np.einsum('i,imn->imn',n[p.nalpha:p.nbf5],J[p.nalpha:p.nbf5,:,:],optimize=True)  # i = [nalpha,nbf5]

    # C^J J
    cj12_ini_nbf5 = view(cj12,ini:p.nbf5,ini:p.nbf5)
    cj12_ini_nbf5[diagind(cj12_ini_nbf5)] .= 0.0
    @tullio F[i,m,n] += cj12[i,j]*J[j,m,n]
    #np.fill_diagonal(cj12[ini:,ini:],0) # Remove diag.
    #F += np.einsum('ij,jmn->imn',cj12,J,optimize=True)                                                # i = [1,nbf5]
    ##F[ini:p.nbf5,:,:] -= np.einsum('ii,imn->imn',cj12[ini:p.nbf5,ini:p.nbf5],J[ini:p.nbf5,:,:],optimize=True) # quita i==j

    # -C^K K
    ck12_ini_nbf5 = view(ck12,ini:p.nbf5,ini:p.nbf5)
    ck12_ini_nbf5[diagind(ck12_ini_nbf5)] .= 0.0
    @tullio F[i,m,n] += -ck12[i,j]*K[j,m,n]
    #np.fill_diagonal(ck12[ini:,ini:],0) # Remove diag.
    #F -= np.einsum('ij,jmn->imn',ck12,K,optimize=True)                                                # i = [1,nbf5]
    ##F[ini:p.nbf5,:,:] += np.einsum('ii,imn->imn',ck12[ini:p.nbf5,ini:p.nbf5],K[ini:p.nbf5,:,:],optimize=True) # quita i==j

    return F
    end

    function computeLagrange(F,C,p)

    Cnbf5 = view(C,:,1:p.nbf5)
    Cnoptorb = view(C,:,1:p.noptorb)

    @tullio G[m,i] := F[i,m,n]*Cnbf5[n,i]
    #G = np.einsum('imn,ni->mi',F,C[:,0:p.nbf5],optimize=True)

    #Compute Lagrange multipliers
    elag = zeros(p.nbf,p.nbf)
    elag_noptorb_nbf5 = view(elag,1:p.noptorb,1:p.nbf5)
    @tullio elag_noptorb_nbf5[i,j] = Cnoptorb[m,i]*G[m,j]
    #elag[0:p.noptorb,0:p.nbf5] = np.einsum('mi,mj->ij',C[:,0:p.noptorb],G,optimize=True)[0:p.noptorb,0:p.nbf5]
    #println(elag)
    return elag

    end

function computeE_elec(H,C,n,elag,p)
    #EELECTRr
    E = 0
    #println("1: ",E)
    elag_nbf5 = view(elag,1:p.nbf5,1:p.nbf5)
    @tullio E += elag_nbf5[i,i]
    #println("2: ",E)
    #E = E + np.einsum('ii',elag[:p.nbf5,:p.nbf5],optimize=True)
    n_beta = view(n,1:p.nbeta)
    C_beta = view(C,:,1:p.nbeta)
    @tullio E += n_beta[i]*C_beta[m,i]*H[m,n]*C_beta[n,i]
    #println("3: ",E)
    #E = E + np.einsum('i,mi,mn,ni',n[:p.nbeta],C[:,:p.nbeta],H,C[:,:p.nbeta],optimize=True)
    if(!p.HighSpin)
        n_beta_alpha = view(n,p.nbeta+1:p.nalpha)
        C_beta_alpha = view(C,:,p.nbeta+1:p.nalpha)
	@tullio E += n_beta_alpha[i]*C_beta_alpha[m,i]*H[m,n]*C_beta_alpha[n,i]
        #println("4: ",E)
    #    E = E + np.einsum('i,mi,mn,ni',n[p.nbeta:p.nalpha],C[:,p.nbeta:p.nalpha],H,C[:,p.nbeta:p.nalpha],optimize=True)
    elseif(p.HighSpin)
        C_beta_alpha = view(C,:,p.nbeta+1:p.nalpha)
	@tullio E += 0.5*C_beta_alpha[m,i]*H[m,n]*C_beta_alpha[n,i]
    #    E = E + 0.5*np.einsum('mi,mn,ni',C[:,p.nbeta:p.nalpha],H,C[:,p.nbeta:p.nalpha],optimize=True)
    end
    n_alpha_nbf5 = view(n,p.nalpha+1:p.nbf5)
    C_alpha_nbf5 = view(C,:,p.nalpha+1:p.nbf5)
    @tullio E += n_alpha_nbf5[i]*C_alpha_nbf5[m,i]*H[m,n]*C_alpha_nbf5[n,i]
    #println("5: ",E)
    #E = E + np.einsum('i,mi,mn,ni',n[p.nalpha:p.nbf5],C[:,p.nalpha:p.nbf5],H,C[:,p.nalpha:p.nbf5],optimize=True)

    return E

    end

function computeLagrangeConvergency(elag)
    # Convergency
    sumdiff = sum(abs.(elag-Transpose(elag)))
    maxdiff = maximum(abs.(elag-Transpose(elag)))

    return sumdiff,maxdiff

end

function fmiug_scaling(fmiug0,elag,i_ext,nzeros,nbf,noptorb)

    #scaling
    fmiug = zeros(nbf,nbf)
    fmiug_noptorb_noptorb = view(fmiug,1:noptorb,1:noptorb)
    if i_ext == 1 && fmiug0==nothing
	@tullio fmiug[i,j] = (elag[i,j] + elag[j,i])/2
#        fmiug[:noptorb,:noptorb] = ((elag[:noptorb,:noptorb] + elag[:noptorb,:noptorb].T) / 2)
    else
	@tullio fmiug[i,j] = (elag[i,j] - elag[j,i])
#        fmiug[:noptorb,:noptorb] = (elag[:noptorb,:noptorb] - elag[:noptorb,:noptorb].T)
        fmiug = tril(fmiug,-1) + Transpose(tril(fmiug,-1))
#        fmiug = np.tril(fmiug,-1) + np.tril(fmiug,-1).T
        for k in 0:nzeros+9
            fmiug[10.0^(9-k) .< abs.(fmiug) .< 10.0^(10-k)] .*= 0.1
        end
#        for k in range(nzeros+9+1):
#            fmiug[(abs(fmiug) > 10**(9-k)) & (abs(fmiug) < 10**(10-k))] *= 0.1
        fmiug[diagind(fmiug)] .= fmiug0            
#        np.fill_diagonal(fmiug[:noptorb,:noptorb],fmiug0[:noptorb])
    end

    return fmiug
    end
    
function fmiug_diis(fk,fmiug,idiis,bdiis,cdiis,maxdiff,p)

    idiis = idiis + 1
    fk[idiis,1:p.noptorb,1:p.noptorb] = fmiug[1:p.noptorb,1:p.noptorb]
    for m in 1:idiis
        bdiis[m,idiis] = 0
        for i in 1:p.noptorb
            for j in 1:i
                bdiis[m,idiis] += fk[m,i,j]*fk[idiis,j,i]
	    end
	end
        bdiis[idiis,m] = bdiis[m,idiis]
        bdiis[m,idiis+1] = -1
        bdiis[idiis+1,m] = -1
    end
    bdiis[idiis+1,idiis+1] = 0

    if idiis>=p.ndiis
        cdiis = zeros(idiis+1)
        cdiis[1:idiis] .= 0
        cdiis[idiis+1] .= -1
        x = bdiis[1:idiis+1,1:idiis+1]\cdiis[1:idiis+1]

        for i in 1:p.noptorb
            for j in 1:i
                fmiug[i,j] = 0
                for k in 1:idiis+1
                    fmiug[i,j] += x[k]*fk[k,i,j]
		end
                fmiug[j,i] = fmiug[i,j]
             end
	end
    end

    if p.perdiis
        idiis = 0
    end
    #else
    #    idiis = idiis + 1
    #end

    return fk,fmiug,idiis,bdiis
end

end
