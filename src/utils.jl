function compute_Lagrange2(C,n,H,I,b_mnl,cj12,ck12,pa)

    Cnbf5 = view(C,1:pa.nbf,1:pa.nbf5)
    @tullio tmp[m,j] := H[m,nn] * Cnbf5[nn,j]
    @tullio H_mat[i,j] := C[m,i] * tmp[m,j]
    if pa.RI
        @tullio tmp[m,q,l] := Cnbf5[nn,q]*b_mnl[m,nn,l]
        @tullio b_MO[p,q,l] := C[m,p]*tmp[m,q,l]
    else
        @tullio tmp[m,q,s,l] := Cnbf5[nn,q]*I[m,nn,s,l]
        @tullio tmp2[m,q,r,l] := Cnbf5[s,r]*tmp[m,q,s,l]
        @tullio tmp[m,q,r,t] := Cnbf5[l,t]*tmp2[m,q,r,l]
        @tullio I_MO[p,q,r,t] := C[m,p]*tmp[m,q,r,t]

    end

    elag = zeros(pa.nbf,pa.nbf)
    cj12[diagind(cj12)] .= 0
    ck12[diagind(ck12)] .= 0

    n_beta =        view(n,1:pa.nbeta)
    n_alpha =       view(n,pa.nalpha+1:pa.nbf5)
    Hmat_nbf5 =      view(H_mat,1:pa.nbf,1:pa.nbf5)
    grad_nbf5 =      view(elag,1:pa.nbf,1:pa.nbf5)
    grad_nbeta =     view(elag,1:pa.nbf,1:pa.nbeta)
    grad_nalpha =    view(elag,1:pa.nbf,pa.nalpha+1:pa.nbf5)
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

    if pa.RI
            if(pa.MSpin==0)
                # 2ndH/dy_ab
                @tullio grad_nbf5[a,b]  += n[b]*Hmat_nbf5[a,b]

                # dJ_pp/dy_ab
                @tullio grad_nbeta[a,b] += n_beta[b]*b_nbf_beta[a,b,k]*b_nbeta_beta[b,b,k]
                @tullio grad_nalpha[a,b] += n_alpha[b]*b_nbf_alpha[a,b,k]*b_nalpha_alpha[b,b,k]

                # C^J_pq dJ_pq/dy_ab 
                @tullio tmp[b,k] := cj12[b,q]*b_nbf5_nbf5[q,q,k]
                @tullio grad_nbf5[a,b] += b_nbf_nbf5[a,b,k]*tmp[b,k]

                # -C^K_pq dK_pq/dy_ab 
                @tullio grad_nbf5[a,b] += -ck12[b,q]*b_nbf_nbf5[a,q,k]*b_nbf5_nbf5[b,q,k]
            end
        else
            if(pa.MSpin==0)
                # 2ndH/dy_ab
                @tullio grad_nbf5[a,b]  += n[b]*Hmat_nbf5[a,b]

                # dJ_pp/dy_ab
                @tullio grad_nbeta[a,b] += n_beta[b]*I_nb_nb_nb[a,b,b,b]
                @tullio grad_nalpha[a,b] += n_alpha[b]*I_na_na_na[a,b,b,b]

                # C^J_pq dJ_pq/dy_ab
                @tullio grad_nbf5[a,b] += cj12[b,q]*I_nbf5_nbf5_nbf5[a,b,q,q]

                # -C^K_pq dK_pq/dy_ab 
                @tullio grad_nbf5[a,b] += -ck12[b,q]*I_nbf5_nbf5_nbf5[a,q,b,q]
        end
    end

    return elag,H_mat

end

function compute_Lagrange2(C,n,H,I,b_mnl::CuArray,cj12,ck12,pa)

    if(pa.gpu)
        C = CuArray{typeof(b_mnl).parameters[1]}(C)
        H = CuArray{typeof(b_mnl).parameters[1]}(H)
    end

    Cnbf5 = C[1:pa.nbf,1:pa.nbf5]
    H_mat = C' * (H * Cnbf5)
    #@tullio H_mat[i,j] := Cnew[m,i] * H[m,nn] * Cnbf5[nn,j]
    if pa.RI
        @tensor tmp[m,q,l] := Cnbf5[nn,q]*b_mnl[m,nn,l]
        @tensor b_MO[p,q,l] := C[m,p]*tmp[m,q,l]
	CUDA.unsafe_free!(tmp)
    else
        @tullio tmp[m,q,s,l] := Cnbf5[nn,q]*I[m,nn,s,l]
        @tullio tmp2[m,q,r,l] := Cnbf5[s,r]*tmp[m,q,s,l]
        @tullio tmp[m,q,r,t] := Cnbf5[l,t]*tmp2[m,q,r,l]
        @tullio I_MO[p,q,r,t] := C[m,p]*tmp[m,q,r,t]

    end

    cj12[diagind(cj12)] .= 0
    ck12[diagind(ck12)] .= 0

    if pa.gpu
        elag = CUDA.zeros(typeof(b_mnl).parameters[1],pa.nbf,pa.nbf)
	cj12 = CuArray{typeof(b_mnl).parameters[1]}(cj12)
	ck12 = CuArray{typeof(b_mnl).parameters[1]}(ck12)
	n = CuArray{typeof(b_mnl).parameters[1]}(n)
    end

    n_beta =        view(n,1:pa.nbeta)
    n_alpha =       view(n,pa.nalpha+1:pa.nbf5)
    Hmat_nbf5 =      view(H_mat,1:pa.nbf,1:pa.nbf5)
    grad_nbf5 =      view(elag,1:pa.nbf,1:pa.nbf5)
    grad_nbeta =     view(elag,1:pa.nbf,1:pa.nbeta)
    grad_nalpha =    view(elag,1:pa.nbf,pa.nalpha+1:pa.nbf5)
    if !pa.RI
        I_nb_nb_nb =    view(I_MO,1:pa.nbf,1:pa.nbeta,1:pa.nbeta,1:pa.nbeta)
        I_na_na_na =    view(I_MO,1:pa.nbf,pa.nalpha+1:pa.nbf5,pa.nalpha+1:pa.nbf5,pa.nalpha+1:pa.nbf5)
        I_nbf5_nbf5_nbf5 =    view(I_MO,1:pa.nbf,1:pa.nbf5,1:pa.nbf5,1:pa.nbf5)
    end

    if pa.RI
            if(pa.MSpin==0)
                # 2ndH/dy_ab
		elag[1:pa.nbf,1:pa.nbf5]  += Hmat_nbf5 .* n'
                #@tullio grad_nbf5[a,b]  += n[b]*Hmat_nbf5[a,b]

                # dJ_pp/dy_ab
		b_nbeta_beta =   view(b_MO,1:pa.nbeta,1:pa.nbeta,1:pa.nbfaux)
                b_nbf_beta =     view(b_MO,1:pa.nbf,1:pa.nbeta,1:pa.nbfaux)

                @tullio tmp[b,k] := b_nbeta_beta[b,b,k]
                tmp2 = n_beta .* tmp
	        CUDA.unsafe_free!(tmp)
		tmp3 = permutedims(b_nbf_beta,(2,3,1))
		tmp4 = tmp2 .* tmp3
	        CUDA.unsafe_free!(tmp2)
	        CUDA.unsafe_free!(tmp3)
		tmp5 = sum(tmp4, dims=2)
	        CUDA.unsafe_free!(tmp4)
		tmp6 = dropdims(tmp5,dims=2)
	        CUDA.unsafe_free!(tmp5)
		elag[1:pa.nbf,1:pa.nbeta] += tmp6'
	        CUDA.unsafe_free!(tmp6)
                #@tullio grad_nbeta[a,b] += n_beta[b]*b_nbf_beta[a,b,k]*b_nbeta_beta[b,b,k]

                b_nalpha_alpha = view(b_MO,pa.nalpha+1:pa.nbf5,pa.nalpha+1:pa.nbf5,1:pa.nbfaux)
                b_nbf_alpha =    view(b_MO,1:pa.nbf,pa.nalpha+1:pa.nbf5,1:pa.nbfaux)
                @tullio tmp[b,k] := b_nalpha_alpha[b,b,k]
                tmp2 = n_alpha .* tmp
                CUDA.unsafe_free!(tmp)
                tmp3 = permutedims(b_nbf_alpha,(2,3,1))
                tmp4 = tmp2 .* tmp3
                CUDA.unsafe_free!(tmp2)
                CUDA.unsafe_free!(tmp3)
                tmp5 = sum(tmp4, dims=2)
                CUDA.unsafe_free!(tmp4)
                tmp6 = dropdims(tmp5,dims=2)
                CUDA.unsafe_free!(tmp5)
                elag[1:pa.nbf,pa.nalpha+1:pa.nbf5] += tmp6'
                CUDA.unsafe_free!(tmp6)
                #@tullio grad_nalpha[a,b] += n_alpha[b]*b_nbf_alpha[a,b,k]*b_nalpha_alpha[b,b,k]

                b_nbf_nbf5 =     view(b_MO,1:pa.nbf,1:pa.nbf5,1:pa.nbfaux)
                b_nbf5_nbf5 =    view(b_MO,1:pa.nbf5,1:pa.nbf5,1:pa.nbfaux)

                # C^J_pq dJ_pq/dy_ab
                @tullio tmp[q,k] := b_nbf5_nbf5[q,q,k]
                @tensor tmp2[b,k] := cj12[b,q]*tmp[q,k]
	        CUDA.unsafe_free!(tmp)
                @tullio grad_nbf5[a,b] += b_nbf_nbf5[a,b,k]*tmp2[b,k]
	        CUDA.unsafe_free!(tmp2)
                #@tullio grad_nbf5[a,b] += b_nbf_nbf5[a,b,k]*cj12[b,q]*b_nbf5_nbf5[q,q,k]

                # -C^K_pq dK_pq/dy_ab
		tmp = ck12 .* b_nbf5_nbf5
		tmp2 = permutedims(tmp,(2,1,3))
	        CUDA.unsafe_free!(tmp)
		tmp3 = NNlibCUDA.batched_mul(b_nbf_nbf5,tmp2)
	        CUDA.unsafe_free!(tmp2)
	        tmp4 = sum(tmp3, dims=3)
		CUDA.unsafe_free!(tmp3)
		elag[1:pa.nbf,1:pa.nbf5] += -dropdims(tmp4, dims=3)
	        CUDA.unsafe_free!(tmp4)
		#@tullio grad_nbf5[a,b] += -ck12[b,q]*b_nbf_nbf5[a,q,k]*b_nbf5_nbf5[b,q,k]
            end
        else
            if(pa.MSpin==0)
                # 2ndH/dy_ab
                @tullio grad_nbf5[a,b]  += n[b]*Hmat_nbf5[a,b]

                # dJ_pp/dy_ab
                @tullio grad_nbeta[a,b] += n_beta[b]*I_nb_nb_nb[a,b,b,b]
                @tullio grad_nalpha[a,b] += n_alpha[b]*I_na_na_na[a,b,b,b]

                # C^J_pq dJ_pq/dy_ab
                @tullio grad_nbf5[a,b] += cj12[b,q]*I_nbf5_nbf5_nbf5[a,b,q,q]

                # -C^K_pq dK_pq/dy_ab
                @tullio grad_nbf5[a,b] += -ck12[b,q]*I_nbf5_nbf5_nbf5[a,q,b,q]
        end
    end

    elag = Array(elag)
    H_mat = Array(H_mat)
    
    return elag,H_mat

end

function ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)

    if(p.no1==0)
        elag,Hmat = compute_Lagrange2(C,n,H,I,b_mnl,cj12,ck12,p)
	E = computeE_elec(Hmat,n,elag,p)
    else
        J,K = computeJKj(C,I,b_mnl,p)

        if p.MSpin==0
            if p.gpu
                F = computeF_RC(J,K,n,H,cj12,ck12,p)
            else
                F = computeF_RC(J,K,n,H,cj12,ck12,p)
            end
        elseif p.MSpin!=0
            F = computeF_RO(J,K,n,H,cj12,ck12,p)
        end

        elag = computeLagrange(F,C,p)
    end
    E = computeE_elec(H,C,n,elag,p)

    sumdiff,maxdiff = computeLagrangeConvergency(elag)

    return E,elag,sumdiff,maxdiff


end

function computeF_RC(J,K,n,H,cj12,ck12,p)

    # Matriz de Fock Generalizada
    F = zeros(p.nbf5,p.nbf,p.nbf)

    ini = 0
    if(p.no1>1)
        ini = p.no1
    end

    # nH
    @tullio F[i,m,s] += n[i]*H[m,s]

    # nJ
    F_ini_beta = view(F,ini+1:p.nbeta,:,:)
    n_ini_beta = view(n,ini+1:p.nbeta)
    J_ini_beta = view(J,ini+1:p.nbeta,:,:)
    @tullio F_ini_beta[i,m,n] += n_ini_beta[i]*J_ini_beta[i,m,n]
    F_alpha_nbf5 = view(F,p.nalpha+1:p.nbf5,:,:)
    n_alpha_nbf5 = view(n,p.nalpha+1:p.nbf5)
    J_alpha_nbf5 = view(J,p.nalpha+1:p.nbf5,:,:)
    @tullio F_alpha_nbf5[i,m,n] += n_alpha_nbf5[i]*J_alpha_nbf5[i,m,n]

    # C^J J
    cj12_ini_nbf5 = view(cj12,ini+1:p.nbf5,ini+1:p.nbf5)
    cj12_ini_nbf5[diagind(cj12_ini_nbf5)] .= 0.0
    @tullio F[i,m,n] += cj12[i,j]*J[j,m,n]

    # -C^K K
    ck12_ini_nbf5 = view(ck12,ini+1:p.nbf5,ini+1:p.nbf5)
    ck12_ini_nbf5[diagind(ck12_ini_nbf5)] .= 0.0
    @tullio F[i,m,n] += -ck12[i,j]*K[j,m,n]

    return F

end

function computeF_RC(J::CuArray,K,n,H,cj12,ck12,p)

    ini = 0
    if(p.no1>1)
        ini = p.no1
    end

    # Matriz de Fock Generalizada
    n = CuArray{typeof(J).parameters[1]}(n)
    H = CuArray{typeof(J).parameters[1]}(H)

    # nH
    @tullio F[i,m,s] := n[i]*H[m,s]

    # -C^K K
    ck12_ini_nbf5 = view(ck12,ini+1:p.nbf5,ini+1:p.nbf5)
    ck12_ini_nbf5[diagind(ck12_ini_nbf5)] .= 0.0
    ck12_gpu = CuArray{typeof(J).parameters[1]}(ck12)
    @tensor F[i,m,n] += -ck12_gpu[i,j]*K[j,m,n]
    CUDA.unsafe_free!(ck12_gpu)
    CUDA.unsafe_free!(K)

    # C^J J
    cj12_ini_nbf5 = view(cj12,ini+1:p.nbf5,ini+1:p.nbf5)
    cj12_ini_nbf5[diagind(cj12_ini_nbf5)] .= 0.0
    cj12_gpu = CuArray{typeof(J).parameters[1]}(cj12)
    @tensor F[i,m,n] += cj12_gpu[i,j]*J[j,m,n]
    CUDA.unsafe_free!(cj12_gpu)

    # nJ
    n_ini_beta = n[ini+1:p.nbeta]
    J_ini_beta = J[ini+1:p.nbeta,1:p.nbf,1:p.nbf]
    F[ini+1:p.nbeta,1:p.nbf,1:p.nbf] += n_ini_beta .* J_ini_beta
    CUDA.unsafe_free!(n_ini_beta)
    CUDA.unsafe_free!(J_ini_beta)

    n_alpha_nbf5 = n[p.nalpha+1:p.nbf5]
    J_alpha_nbf5 = J[p.nalpha+1:p.nbf5,1:p.nbf,1:p.nbf]
    F[p.nalpha+1:p.nbf5,1:p.nbf,1:p.nbf] += n_alpha_nbf5 .* J_alpha_nbf5
    CUDA.unsafe_free!(n_alpha_nbf5)
    CUDA.unsafe_free!(J_alpha_nbf5)
    CUDA.unsafe_free!(J)

    return Array(F)

end

function computeLagrange(F,C,p)

    Cnbf5 = view(C,:,1:p.nbf5)
    Cnoptorb = view(C,:,1:p.noptorb)

    @tullio G[m,i] := F[i,m,n]*Cnbf5[n,i]

    #Compute Lagrange multipliers
    elag = zeros(p.nbf,p.nbf)
    elag_noptorb_nbf5 = view(elag,1:p.noptorb,1:p.nbf5)
    @tullio elag_noptorb_nbf5[i,j] = Cnoptorb[m,i]*G[m,j]

    return elag

    end

function computeF_RO(J,K,n,H,cj12,ck12,p)

    # Matriz de Fock Generalizada
    F = zeros(p.nbf5,p.nbf,p.nbf)

    ini = 0
    if p.no1>1
        ini = p.no1
    end

    F_ini_beta = view(F,ini+1:p.nbeta,:,:)
    F_beta = view(F,1:p.nbeta,:,:)
    F_alpha = view(F,p.nbeta+1:p.nalpha,:,:)
    F_nbf5 = view(F,p.nalpha+1:p.nbf5,:,:)
    n_ini_beta = view(n,ini+1:p.nbeta)
    n_beta = view(n,1:p.nbeta)
    n_alpha = view(n,p.nbeta+1:p.nalpha)
    n_nbf5 = view(n,p.nalpha+1:p.nbf5)
    J_ini_beta = view(J,ini+1:p.nbeta,:,:)
    J_beta = view(J,1:p.nbeta,:,:)
    J_alpha = view(J,p.nbeta+1:p.nalpha,:,:)
    J_nbf5 = view(J,p.nalpha+1:p.nbf5,:,:)
    K_beta = view(K,1:p.nbeta,:,:)
    K_alpha = view(K,p.nbeta+1:p.nalpha,:,:)
    K_nbf5 = view(K,p.nalpha+1:p.nbf5,:,:)

    # nH
    @tullio F_beta[i,m,n] += n_beta[i]*H[m,n]
    @tullio F_alpha[i,m,n] += 0.5*H[m,n]
    @tullio F_nbf5[i,m,n] += n_nbf5[i]*H[m,n]

    # nJ
    @tullio avx=false F_ini_beta[i,m,n] += n_ini_beta[i]*J_ini_beta[i,m,n]
    @tullio F_nbf5[i,m,n] += n_nbf5[i]*J_nbf5[i,m,n]

    # C^J J
    cj12_ini_nbf5 = view(cj12,ini+1:p.nbf5,ini+1:p.nbf5)
    cj12_ini_nbf5[diagind(cj12_ini_nbf5)] .= 0.0
    cj12_beta_beta = view(cj12,1:p.nbeta,1:p.nbeta)
    cj12_beta_nbf5 = view(cj12,1:p.nbeta,p.nalpha+1:p.nbf5)
    cj12_nbf5_beta = view(cj12,p.nalpha+1:p.nbf5,1:p.nbeta)
    cj12_nbf5_nbf5 = view(cj12,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
    @tullio F_beta[i,m,n] += cj12_beta_beta[i,j]*J_beta[j,m,n]
    @tullio F_beta[i,m,n] += cj12_beta_nbf5[i,j]*J_nbf5[j,m,n]
    @tullio F_nbf5[i,m,n] += cj12_nbf5_beta[i,j]*J_beta[j,m,n]
    @tullio F_nbf5[i,m,n] += cj12_nbf5_nbf5[i,j]*J_nbf5[j,m,n]

    # -C^K K
    ck12_ini_nbf5 = view(ck12,ini+1:p.nbf5,ini+1:p.nbf5)
    ck12_ini_nbf5[diagind(ck12_ini_nbf5)] .= 0.0
    ck12_beta_beta = view(ck12,1:p.nbeta,1:p.nbeta)
    ck12_beta_nbf5 = view(ck12,1:p.nbeta,p.nalpha+1:p.nbf5)
    ck12_nbf5_beta = view(ck12,p.nalpha+1:p.nbf5,1:p.nbeta)
    ck12_nbf5_nbf5 = view(ck12,p.nalpha+1:p.nbf5,p.nalpha+1:p.nbf5)
    @tullio F_beta[i,m,n] += -ck12_beta_beta[i,j]*K_beta[j,m,n]
    @tullio F_beta[i,m,n] += -ck12_beta_nbf5[i,j]*K_nbf5[j,m,n]
    @tullio F_nbf5[i,m,n] += -ck12_nbf5_beta[i,j]*K_beta[j,m,n]
    @tullio F_nbf5[i,m,n] += -ck12_nbf5_nbf5[i,j]*K_nbf5[j,m,n]

    # SUMij
    @tullio F_beta[i,m,n] += n_beta[i]*J_alpha[j,m,n]
    @tullio F_beta[i,m,n] += -0.5*n_beta[i]*K_alpha[j,m,n]
    @tullio avx=false F_alpha[i,m,n] += 0.5*J_alpha[j,m,n]
    @tullio avx=false F_alpha[i,m,n] += -0.5*K_alpha[j,m,n]
    F[p.nbeta+1:p.nalpha,:,:] -= 0.5*(J[p.nbeta+1:p.nalpha,:,:]-K[p.nbeta+1:p.nalpha,:,:]) #Remove diag.
    @tullio F_nbf5[i,m,n] += n_nbf5[i]*J_alpha[j,m,n]
    @tullio F_nbf5[i,m,n] += -0.5*n_nbf5[i]*K_alpha[j,m,n]

    # PRODWROij
    @tullio avx=false F_alpha[i,m,n] += n_beta[j]*J_beta[j,m,n]
    @tullio avx=false F_alpha[i,m,n] += -0.5*n_beta[j]*K_beta[j,m,n]
    @tullio avx=false F_alpha[i,m,n] += n_nbf5[j]*J_nbf5[j,m,n]
    @tullio avx=false F_alpha[i,m,n] += -0.5*n_nbf5[j]*K_nbf5[j,m,n]

    return F
end

function computeE_elec(H,C,n,elag,p)
    #EELECTRr
    E = 0
    elag_nbf5 = view(elag,1:p.nbf5,1:p.nbf5)
    @tullio E += elag_nbf5[i,i]
    n_beta = view(n,1:p.nbeta)
    C_beta = view(C,:,1:p.nbeta)
    @tullio E += n_beta[i]*C_beta[m,i]*H[m,n]*C_beta[n,i]
    if(!p.HighSpin)
        n_beta_alpha = view(n,p.nbeta+1:p.nalpha)
        C_beta_alpha = view(C,:,p.nbeta+1:p.nalpha)
	@tullio avx=false E += n_beta_alpha[i]*C_beta_alpha[m,i]*H[m,n]*C_beta_alpha[n,i]
    elseif(p.HighSpin)
        C_beta_alpha = view(C,:,p.nbeta+1:p.nalpha)
	@tullio avx=false E += 0.5*C_beta_alpha[m,i]*H[m,n]*C_beta_alpha[n,i]
    end
    n_alpha_nbf5 = view(n,p.nalpha+1:p.nbf5)
    C_alpha_nbf5 = view(C,:,p.nalpha+1:p.nbf5)
    @tullio E += n_alpha_nbf5[i]*C_alpha_nbf5[m,i]*H[m,n]*C_alpha_nbf5[n,i]

    return E

    end

function computeE_elec(Hmat,n,elag,p)
    #EELECTRr
    E = 0
    elag_nbf5 = view(elag,1:p.nbf5,1:p.nbf5)
    @tullio E += elag_nbf5[i,i]
    n_nbf5 = view(n,1:p.nbf5)
    H_nbf5 = view(Hmat,1:p.nbf5,1:p.nbf5)
    @tullio E += n_nbf5[i]*H_nbf5[i,i]

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
    if i_ext == 1 && isnothing(fmiug0)
	@tullio fmiug[i,j] = (elag[i,j] + elag[j,i])/2
    else
	@tullio fmiug[i,j] = (elag[i,j] - elag[j,i])
        fmiug = tril(fmiug,-1) + Transpose(tril(fmiug,-1))
        for k in 0:nzeros+9
            fmiug[10.0^(9-k) .< abs.(fmiug) .< 10.0^(10-k)] .*= 0.1
        end
        fmiug[diagind(fmiug)] .= fmiug0            
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

    return fk,fmiug,idiis,bdiis
end

function compute_E_nuc(bset,p)

    E_nuc = 0.0
    for i in 1:p.natoms
        for j in i+1:p.natoms
            E_nuc += bset.atoms[i].Z*bset.atoms[j].Z/(norm(bset.atoms[i].xyz-bset.atoms[j].xyz)*1.88973)
        end
    end

    return E_nuc

end

function Z_to_symbol(Z)
    dict = Dict(1 => "H", 2 => "He", 3 => "Li",
        4 => "Be", 5 => "B", 6 => "C", 7 => "N", 8 => "O", 9 => "F", 10 => "Ne",
        11 => "Na", 12 => "Mg", 13 => "Al", 14 => "Si", 15 => "P", 16 => "S", 17 => "Cl", 18 => "Ar",
        19 => "K", 20 => "Ca", 21 => "Sc", 22 => "Ti", 23 => "V", 24 => "Cr", 25 => "Mn", 26 => "Fe", 27 => "Co", 28 => "Ni", 29 => "Cu", 30 => "Zn", 31 => "Ga", 32 => "Ge", 33 => "As", 34 => "Se", 35 => "Br", 36 => "Kr",
        37 => "Rb", 38 => "Sr", 39 => "Y", 40 => "Zr", 41 => "Nb", 42 => "Mo", 43 => "Tc", 44 => "Ru", 45 => "Rh", 46 => "Pd", 47 => "Ag", 48 => "Cd", 49 => "In", 50 => "Sn", 51 => "Sb", 52 => "Te", 53 => "I", 54 => "Xe")

    return dict[Z]
end

function check_ortho(C,S,p)

    # Revisa ortonormalidad
    orthonormality = true
    CTSC = C'*S*C
    ortho_deviation = abs.(CTSC - I)
    if (maximum(ortho_deviation) > 10^-6)
        orthonormality = false
    end
    if !orthonormality
        @printf("Orthonormality violations %i, Maximum Violation %f\n",sum(ortho_deviation .> 10^-6),maximum(ortho_deviation))
        println("Trying to orthonormalize")
        C = orthonormalize(C,S,p)
        C = check_ortho(C,S,p)
    else
        println("No violations of the orthonormality")
    end
    for j in 1:p.nbf
        #Obtiene el índice del coeficiente con mayor valor absoluto del MO
        idxmaxabsval = 1
        for i in 1:p.nbf
	    if(abs(C[i,j])>abs(C[idxmaxabsval,j]))
                 idxmaxabsval = i
            end
        end
        # Ajusta el signo del MO
        C[1:p.nbf,j] = sign(C[idxmaxabsval,j])*C[1:p.nbf,j]
    end
        
    return C

end

function orthonormalize(C,S,p)

    evals,evecs = eigen(S)
    for i in 1:size(evals)[1]
        if (evals[i]<0.0)
            evals[i] = 0.0
        else
            evals[i] = 1/sqrt(evals[i])
        end
    end

    @tullio S_12[i,j] := evecs[i,j]*evals[j]

    @tullio Cnew[i,j] := S[i,k]*C[k,j] #np.einsum('ik,kj->ij',S,C,optimize=True)

    @tullio Cnew2[i,j] := S_12[k,i]*Cnew[k,j]  #np.einsum('ki,kj->ij',S_12,Cnew)

    for i in 1:size(Cnew2)[2]
	Cnew2[1:p.nbf,i] = Cnew2[1:p.nbf,i]/norm(Cnew2[1:p.nbf,i])
	for j in i+1:size(Cnew2)[1]
	    val = -sum(Cnew2[1:p.nbf,i].*Cnew2[1:p.nbf,j])
            Cnew2[1:p.nbf,j] = Cnew2[1:p.nbf,j] + val*Cnew2[1:p.nbf,i]
        end
    end

    @tullio C[i,j] = S_12[i,k]*Cnew2[k,j]

    return C

end

function rotate_orbital(y,C,p)

    ynew = zeros(p.nbf,p.nbf)

    n = 1
    for i in 1:p.nbf5
        for j in i+1:p.nbf
            ynew[i,j] =  y[n]
            ynew[j,i] = -y[n]
            n += 1
	end
    end

    U = exp(ynew)
    U = real.(U)
    @tullio Cnew[m,p] := C[m,r]*U[r,p]

    return Cnew

end

function n_to_gammas_softmax(n,p)

    gamma = zeros(p.nv)

    for i in 1:p.ndoc

        ll = p.no1 + p.ndns + p.ncwo*(p.ndoc - i) + 1
        ul = p.no1 + p.ndns + p.ncwo*(p.ndoc - i + 1)
        llg = ll - p.ndns + p.ndoc - p.no1
        ulg = ul - p.ndns + p.ndoc - p.no1

	ns = n[ll:ul]

	A = zeros(p.ncwo,p.ncwo)
	b = zeros(p.ncwo)

        for j in 1:p.ncwo
	    A[j,:] .= ns[j]
	    A[j,j] = ns[j]-1
	    b[j] = -ns[j]
	end

	x = log.(A\b)

	gamma[llg:ulg] .= x

    end

    return gamma

end

function n_to_gammas_trigonometric(n,p)
    gamma = zeros(p.nv)
    for i in 1:p.ndoc
	idx = p.no1 + i
	gamma[i] = acos(sqrt(2.0*n[idx]-1.0))
	prefactor = max(1-n[idx],1e-14)
        for j in 1:p.ncwo-1
            jg = p.ndoc+(i-1)*(p.ncwo-1)+j
            ig = p.no1 + p.ndns + p.ncwo*(p.ndoc - i) + j
	    gamma[jg] = asin(sqrt(n[ig]/prefactor))
	    prefactor = prefactor * (cos(gamma[jg]))^2
        end
    end
    return gamma
end
