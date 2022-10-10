function nofmp2(n,C,H,I,b_mnl,E_nuc,p,nofmp2strategy,tol)

    println(" NOF-MP2")
    println("=========")

    occ = view(n,p.no1+1:p.nbf5)
    vec = view(C,:,p.no1+1:p.nbf)
    D = computeD_HF(C,I,b_mnl,p)
    if p.MSpin==0
        if p.nsoc>0
            Dalpha = computeDalpha_HF(C,I,b_mnl,p)
            D = D + 0.5*Dalpha
        end
        J,K = computeJK_HF(D,I,p)
        F = H + 2*J - K
	@tullio DH[i,j] := D[i,k]*H[k,j]
	@tullio DF[i,j] := D[i,k]*F[k,j]
	DHDF = DH + DF
	@tullio EHFL := DHDF[i,i]
        #EHFL = np.trace(np.matmul(D,H)+np.matmul(D,F))
    elseif p.MSpin!=0
        D = computeD_HF(C,I,b_mnl,p)
        Dalpha = computeDalpha_HF(C,I,b_mnl,p)
        J,K = computeJK_HF(D,I,b_mnl,p)
        F = 2*J - K
	DDa1 = D + 0.5*Dalpha
	@tullio DDa1H[i,j] := DDa1[i,k]*H[k,j]
	DDa2 = D + Dalpha
	@tullio DDa2F[i,j] := DDa2[i,k]*F[k,j]
	DDa12HF = DDa1H + DDa2F
	@tullio EHFL := 2*DDa12HF[i,i]
        #EHFL = 2*np.trace(np.matmul(D+0.5*Dalpha,H))+np.trace(np.matmul(D+Dalpha,F))
        F = H + F
        if p.nsoc>1
            J,K = computeJK_HF(0.5*Dalpha,I,b_mnl,p)
            Falpha = J - K
	    @tullio DaFa[i,j] := 0.5*Dalpha[i,k]*Falpha[k,j]
	    @tullio EHFL += 2*DaFa[i,i]
            #EHFL = EHFL + 2*np.trace(np.matmul(0.5*Dalpha,Falpha))
            F = F + Falpha
        end
    end

    @tullio F_MO[i,j] := vec[k,i]*F[k,l]*vec[l,j]

    F_MO_act = view(F_MO,1:p.nbf-p.no1,1:p.nbf-p.no1)
    @tullio eig[i] := F_MO_act[i,i]
    #eig = np.einsum("ii->i",F_MO[:p.nbf-p.no1,:p.nbf-p.no1])

    iajb = compute_iajb(C,I,p)
    #iajb = compute_iajb(C,I,b_mnl,p)
    
    FI1 = ones(p.nbf-p.no1)
    FI2 = ones(p.nbf-p.no1)

    FI1[1:p.nbf5-p.no1] = 1.0 .- (1.0 .- abs.(1.0 .-2*occ[1:p.nbf5-p.no1])).^2

    FI2[p.nalpha-p.no1+1:p.nbf5-p.no1] = abs.(1.0 .-2*occ[p.nalpha-p.no1+1:p.nbf5-p.no1]).^2

    Tijab = CalTijab(iajb,F_MO,eig,FI1,FI2,p,nofmp2strategy,tol)
    ECd = 0

    for k in 1:p.nvir
        for l in 1:p.nvir
            for i in 1:p.ndoc
                for j in 1:p.ndoc
                    Xijkl = iajb[j,k,i,l]
                    ijkl = i+(j-1)*p.ndns+(k-1)*p.ndns*p.ndns+(l-1)*p.ndns*p.ndns*p.nvir
                    ijlk = i+(j-1)*p.ndns+(l-1)*p.ndns*p.ndns+(k-1)*p.ndns*p.ndns*p.nvir
                    ECd = ECd + Xijkl*(2*Tijab[ijkl]-Tijab[ijlk])
                end
                for j in p.ndoc+1:p.ndns
                    Xijkl = iajb[j,k,i,l]
                    ijkl = i+(j-1)*p.ndns+(k-1)*p.ndns*p.ndns+(l-1)*p.ndns*p.ndns*p.nvir
                    ijlk = i+(j-1)*p.ndns+(l-1)*p.ndns*p.ndns+(k-1)*p.ndns*p.ndns*p.nvir
                    ECd = ECd + Xijkl*(Tijab[ijkl]-0.5*Tijab[ijlk])
                end
            end
            for i in p.ndoc+1:p.ndns
                for j in 1:p.ndoc
                    Xijkl = iajb[j,k,i,l]
                    ijkl = i+(j-1)*p.ndns+(k-1)*p.ndns*p.ndns+(l-1)*p.ndns*p.ndns*p.nvir
                    ijlk = i+(j-1)*p.ndns+(l-1)*p.ndns*p.ndns+(k-1)*p.ndns*p.ndns*p.nvir
                    ECd = ECd + Xijkl*(Tijab[ijkl]-0.5*Tijab[ijlk])
                end
                for j in p.ndoc+1:p.ndns
                    Xijkl = iajb[j,k,i,l]
                    ijkl = i+(j-1)*p.ndns+(k-1)*p.ndns*p.ndns+(l-1)*p.ndns*p.ndns*p.nvir
                    ijlk = i+(j-1)*p.ndns+(l-1)*p.ndns*p.ndns+(k-1)*p.ndns*p.ndns*p.nvir
                    if(j!=i)
                        ECd = ECd + Xijkl*(Tijab[ijkl]-0.5*Tijab[ijlk])/2
                    end
                end
            end
        end
    end

    fi = 2*n.*(1 .-n)

    @tullio CK12nd[i,j] := fi[i]*fi[j]
    #CK12nd = np.outer(fi,fi)

    beta = sqrt.((1 .-abs.(1 .-2*n)).*n)
    #beta = np.sqrt((1-abs(1-2*n))*n)

    for l in 1:p.ndoc
        ll = p.no1 + p.ndns + p.ncwo*(p.ndoc - l) + 1
        ul = p.no1 + p.ndns + p.ncwo*(p.ndoc - l + 1)
        CK12nd[p.no1+l,ll:ul] .= beta[p.no1+l]*beta[ll:ul]
        CK12nd[ll:ul,p.no1+l] .= beta[ll:ul]*beta[p.no1+l]
        CK12nd[ll:ul,ll:ul] .= - beta[ll:ul] .* transpose(beta[ll:ul])
    end

    #C^K KMO
    J_MO,K_MO,H_core = computeJKH_MO(C,H,I,b_mnl,p)

    ECndHF = 0
    ECndl = 0
    if (p.MSpin==0)
        if p.nbeta != p.nalpha
            CK12nd_alpha_beta = view(CK12nd,p.nbeta:p.nalpha,p.nbeta:p.nalpha)
            K_MO_alpha_beta = view(K_MO,p.nbeta:p.nalpha,p.nbeta:p.nalpha)
            @tullio avx=false ECndHF += -CK12nd_alpha_beta[i,j]*K_MO_alpha_beta[j,i]
        end
    #   ECndHF = - np.einsum('ii,ii->',CK12nd[p.nbeta:p.nalpha,p.nbeta:p.nalpha],K_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha]) # sum_ij
        @tullio avx=false ECndl += -CK12nd[i,j]*K_MO[j,i]
        @tullio avx=false ECndl += CK12nd[i,i]*K_MO[i,i]
    #   ECndl -= np.einsum('ij,ji->',CK12nd,K_MO) # sum_ij
    #   ECndl += np.einsum('ii,ii->',CK12nd,K_MO) # Quita i=j
    elseif (p.MSpin!=0)
        CK12nd_no1_beta_no1_nbeta = view(CK12nd_no1_beta,p.no1:p.nbeta,p.no1:p.nbeta)
        CK12nd_no1_beta_alpha_beta = view(CK12nd_no1_beta,p.no1:p.nbeta,p.nalpha:p.nbeta)
        CK12nd_alpha_beta_no1_beta = view(CK12nd_no1_beta,p.nalpha:p.nbeta,p.no1:p.nbeta)
        CK12nd_alpha_beta_alpha_beta = view(CK12nd_alpha_beta,p.nbeta:p.nalpha,p.nbeta:p.nalpha)
        K_MO_no1_beta_no1_nbeta = view(CK12nd_no1_beta,p.no1:p.nbeta,p.no1:p.nbeta)
        K_MO_no1_beta_alpha_beta = view(CK12nd_no1_beta,p.no1:p.nbeta,p.nalpha:p.nbeta)
        K_MO_alpha_beta_no1_beta = view(CK12nd_no1_beta,p.nalpha:p.nbeta,p.no1:p.nbeta)
        K_MO_alpha_beta_alpha_beta = view(CK12nd_alpha_beta,p.nbeta:p.nalpha,p.nbeta:p.nalpha)
        @tullio avx=false ECndl += -CK12nd_no1_beta_no1_nbeta[i,j]*K_MO_no1_beta_no1_nbeta[j,i]
        @tullio avx=false ECndl += -CK12nd_no1_beta_alpha_beta[i,j]*K_MO_alpha_beta_no1_nbeta[j,i]
        @tullio avx=false ECndl += -CK12nd_alpha_beta_no1_nbeta[i,j]*K_MO_no1_beta_alpha_beta[j,i]
        @tullio avx=false ECndl += -CK12nd_alpha_beta_alpha_beta[i,j]*K_MO_alpha_beta_alpha_beta[j,i]
        @tullio avx=false ECndl += CK12nd_no1_beta_no1_nbeta[i,i]*K_MO_no1_beta_no1_nbeta[i,i]
        @tullio avx=false ECndl += CK12nd_alpha_beta_alpha_beta[i,i]*K_MO_alpha_beta_alpha_beta[i,i]
    #   ECndl -= np.einsum('ij,ji->',CK12nd[p.no1:p.nbeta,p.no1:p.nbeta],K_MO[p.no1:p.nbeta,p.no1:p.nbeta]) # sum_ij
    #   ECndl -= np.einsum('ij,ji->',CK12nd[p.no1:p.nbeta,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,p.no1:p.nbeta]) # sum_ij
    #   ECndl -= np.einsum('ij,ji->',CK12nd[p.nalpha:p.nbf5,p.no1:p.nbeta],K_MO[p.no1:p.nbeta,p.nalpha:p.nbf5]) # sum_ij
    #   ECndl -= np.einsum('ij,ji->',CK12nd[p.nalpha:p.nbf5,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5]) # sum_ij
    #   ECndl += np.einsum('ii,ii->',CK12nd[p.no1:p.nbeta,p.no1:p.nbeta],K_MO[p.no1:p.nbeta,p.no1:p.nbeta]) # Quita i=j
    #   ECndl += np.einsum('ii,ii->',CK12nd[p.nalpha:p.nbf5,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5]) # Quita i=j

    end
 
    @printf("      Ehfc      = %15.7f\n",EHFL+E_nuc+ECndHF)
    @printf("\n")
    @printf("      ECd       = %15.7f\n",ECd)
    @printf("      ECnd      = %15.7f\n",ECndl)
    @printf("      Ecorre    = %15.7f\n",ECd+ECndl)
    @printf("      E(NOFMP2) = %15.7f\n",EHFL+ECd+ECndl+E_nuc+ECndHF)
    @printf("\n")

end

function CalTijab(iajb,F_MO,eig,FI1,FI2,p,nofmp2strategy,tol)

    println("Starting CalTijab")

    B = build_B(iajb,FI1,FI2,p.ndoc,p.ndns,p.nvir,p.ncwo)
    println("....B vector Computed")

    if(nofmp2strategy=="numerical")
        println("....Numerical Strategy for Tijab")
        Tijab = Tijab_guess(iajb,eig,p.ndoc,p.ndns,p.nvir)
        println("........Tijab Guess Computed")
        Tijab = vec(Tijab)
        R_norm = build_R(Tijab,B,F_MO,FI1,FI2,p.no1,p.ndoc,p.ndns,p.nvir,p.ncwo,p.nbf)
        @printf("............norm of the residual vector of Tijab Guess %4.1e\n",R_norm)
        res = optimize(Tijab->build_R(Tijab,B,F_MO,FI1,FI2,p.no1,p.ndoc,p.ndns,p.nvir,p.ncwo,p.nbf),Tijab,LBFGS())
        Tijab = Optim.minimizer(res)
    elseif(nofmp2strategy=="analytical")
        println("....Analytical Strategy for Tijab")
        B = transpose(B)
        A = build_A(F_MO,FI1,FI2,p.no1,p.ndoc,p.ndns,p.nvir,p.ncwo,p.nbf,tol)
        println("........A matrix built")
        @printf("............It has %d/%d elements with Tol = %4.1e\n",nnz(A),p.nvir^4*p.ndoc^4,tol)
	Tijab = B/A
        B = transpose(B)
    end

    println("........Tijab found")

    R_norm = build_R(Tijab,B,F_MO,FI1,FI2,p.no1,p.ndoc,p.ndns,p.nvir,p.ncwo,p.nbf)
    @printf("............norm of the residual vector of Tijab solution %4.1e\n",R_norm)
    println("")

    return Tijab

end

function build_B(iajb,FI1,FI2,ndoc,ndns,nvir,ncwo)
    B = zeros(ndns^2*nvir^2)
    @Threads.threads  for i in 1:ndns
        lmin_i = ndns+ncwo*(ndns-i)+1
        lmax_i = ndns+ncwo*(ndns-i)+ncwo
        for j in 1:ndns
            if i==j
                for k in 1:nvir
		    ik = i + (k-1)*ndns
                    kn = k + ndns
                    for l in 1:nvir
                        ln = l + ndns
                        if lmin_i <= kn && kn <= lmax_i && lmin_i <= ln && ln <= lmax_i
                            Ciikl = FI1[kn]*FI1[ln]*FI1[i]*FI1[i]
                        else
                            Ciikl = FI2[kn]*FI2[ln]*FI2[i]*FI2[i]
			end
			iikl =  i + (i-1)*ndns + (k-1)*ndns*ndns + (l-1)*ndns*ndns*nvir
                        B[iikl] = - Ciikl*iajb[i,k,i,l]
		    end
		end
            else
                for k in 1:nvir
		    ik = i + (k-1)*ndns
                    kn = k + ndns
                    for l in 1:nvir
                        ln = l + ndns
			ijkl =  i + (j-1)*ndns + (k-1)*ndns*ndns + (l-1)*ndns*ndns*nvir
                        Cijkl = FI2[kn]*FI2[ln]*FI2[i]*FI2[j]
                        B[ijkl] = - Cijkl*iajb[j,k,i,l]
		    end
		end
	    end
	end
    end
    return B

end

function Tijab_guess(iajb,eig,ndoc,ndns,nvir)
    Tijab = zeros(nvir^2*ndns^2)
    @Threads.threads for ia in 1:nvir
        for i in 1:ndns
            for ib in 1:nvir
                for j in 1:ndns
                    ijab = i + (j-1)*ndns + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir
                    Eijab = eig[ib+ndns] + eig[ia+ndns] - eig[j] - eig[i]
                    Tijab[ijab] = - iajb[j,ia,i,ib]/Eijab
		end
	    end
	end
    end
    return Tijab

end

function build_R(T,B,F_MO,FI1,FI2,no1,ndoc,ndns,nvir,ncwo,nbf)
    npair = zeros(nvir)
    for i in 1:ndoc
        ll = ncwo*(ndoc - i) + 1
        ul = ncwo*(ndoc - i + 1)
        npair[ll:ul] .= i
    end

    Bp = zeros(ndns^2*nvir^2)
    @Threads.threads for ib in 1:nvir
        for ia in 1:nvir
            for j in 1:ndns
                for i in 1:ndns
                    jab =     (j-1)*ndns + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir
                    iab = i              + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir
                    ijb = i + (j-1)*ndns                    + (ib-1)*ndns*ndns*nvir
                    ija = i + (j-1)*ndns + (ia-1)*ndns*ndns
                    ijab= i + (j-1)*ndns + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir

                    Bp[ijab] += (F_MO[ia+ndns,ia+ndns] + F_MO[ib+ndns,ib+ndns] - F_MO[i,i] - F_MO[j,j])*T[i+jab]

                    for k in 1:i-1
			if abs(F_MO[i,k])>1e-10
                            Cki = FI2[k]*FI2[i]
			    Bp[ijab] += (- Cki*F_MO[i,k])*T[k+jab]
		        end
		    end
                    
                    for k in i+1:ndns
			if abs(F_MO[i,k])>1e-10
                            Cki = FI2[k]*FI2[i]
	    		    Bp[ijab] += (- Cki*F_MO[i,k])*T[k+jab]
		        end
		    end

                    for k in 1:j-1
                        if abs(F_MO[j,k])>1e-10
                            Ckj = FI2[k]*FI2[j]
			    Bp[ijab] += (- Ckj*F_MO[j,k])*T[(k-1)*ndns+iab]
		        end
		    end
                    for k in j+1:ndns
                        if abs(F_MO[j,k])>1e-10
                            Ckj = FI2[k]*FI2[j]
    		            Bp[ijab] += (- Ckj*F_MO[j,k])*T[(k-1)*ndns+iab]
		        end
		    end

                    for k in 1:ia-1
                        if abs(F_MO[ia+ndns,k+ndns])>1e-10
                            if npair[k]==npair[ia]
                                Ckia = FI1[k+ndns]*FI1[ia+ndns]
                            else
                                Ckia = FI2[k+ndns]*FI2[ia+ndns]
		            end
			    Bp[ijab] += (Ckia*F_MO[ia+ndns,k+ndns]) * T[(k-1)*ndns*ndns + ijb]
		        end
		    end
                    for k in ia+1:nvir
                        if abs(F_MO[ia+ndns,k+ndns])>1e-10
                            if npair[k]==npair[ia]
                                Ckia = FI1[k+ndns]*FI1[ia+ndns]
                            else
                                Ckia = FI2[k+ndns]*FI2[ia+ndns]
		            end
			    Bp[ijab] += (Ckia*F_MO[ia+ndns,k+ndns]) * T[(k-1)*ndns*ndns + ijb]
		        end
		    end

                    for k in 1:ib-1
                        if abs(F_MO[ib+ndns,k+ndns])>1e-10
                            if npair[k]==npair[ib]
                                Ckib = FI1[k+ndns]*FI1[ib+ndns]
                            else
                                Ckib = FI2[k+ndns]*FI2[ib+ndns]
		            end
			    Bp[ijab] += (Ckib*F_MO[ib+ndns,k+ndns]) * T[(k-1)*ndns*ndns*nvir + ija]
		        end
		    end
                    for k in ib+1:nvir
                        if abs(F_MO[ib+ndns,k+ndns])>1e-10
                            if npair[k]==npair[ib]
                                Ckib = FI1[k+ndns]*FI1[ib+ndns]
                            else
                                Ckib = FI2[k+ndns]*FI2[ib+ndns]
		            end
			    Bp[ijab] += (Ckib*F_MO[ib+ndns,k+ndns]) * T[(k-1)*ndns*ndns*nvir + ija]
		        end
		    end
		end
	    end
        end
    end

    R = B-Bp

    return norm(R)^2
end



function build_A(F_MO,FI1,FI2,no1,ndoc,ndns,nvir,ncwo,nbf,tol)
    npair = zeros(nvir)
    for i in 1:ndoc
        ll = ncwo*(ndoc - i) + 1
        ul = ncwo*(ndoc - i + 1)
        npair[ll:ul] .= i
    end
    A = spzeros(ndns^2*nvir^2,ndns^2*nvir^2)

    l = Threads.SpinLock()
    @Threads.threads for ib in 1:nvir
        for ia in 1:nvir
            for j in 1:ndns
                for i in 1:ndns
                    jab =     (j-1)*ndns + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir
                    iab = i              + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir
                    ijb = i + (j-1)*ndns                    + (ib-1)*ndns*ndns*nvir
                    ija = i + (j-1)*ndns + (ia-1)*ndns*ndns
                    ijab= i + (j-1)*ndns + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir

		    Threads.lock(l)
                    A[ijab,i + jab] = (F_MO[ia+ndns,ia+ndns] + F_MO[ib+ndns,ib+ndns] - F_MO[i,i] - F_MO[j,j])
		    Threads.unlock(l)

                    for k in 1:i-1
                        if abs(F_MO[i,k])>tol
                            Cki = FI2[k]*FI2[i]
		            Threads.lock(l)
                            A[ijab,k + jab]=(- Cki*F_MO[i,k])
                            A[k+jab,ijab]=(- Cki*F_MO[i,k])
		            Threads.unlock(l)
			end
		    end
                    for k in 1:j-1
                        if abs(F_MO[j,k])>tol
                            Ckj = FI2[k]*FI2[j]
		            Threads.lock(l)
			    A[ijab,(k-1)*ndns+iab]=(- Ckj*F_MO[j,k])
			    A[(k-1)*ndns+iab,ijab]=(- Ckj*F_MO[j,k])
		            Threads.unlock(l)
			end
		    end

                    for k in 1:ia-1
                        if abs(F_MO[ia+ndns,k+ndns])>tol
                            if npair[k]==npair[ia]
                                Ckia = FI1[k+ndns]*FI1[ia+ndns]
                            else
                                Ckia = FI2[k+ndns]*FI2[ia+ndns]
			    end
		            Threads.lock(l)
			    A[ijab,(k-1)*ndns*ndns + ijb]=(Ckia*F_MO[ia+ndns,k+ndns])
			    A[(k-1)*ndns*ndns + ijb,ijab]=(Ckia*F_MO[ia+ndns,k+ndns])
		            Threads.unlock(l)
			end
		    end

                    for k in 1:ib-1
                        if abs(F_MO[ib+ndns,k+ndns])>tol
                            if npair[k]==npair[ib]
                                Ckib = FI1[k+ndns]*FI1[ib+ndns]
                            else
                                Ckib = FI2[k+ndns]*FI2[ib+ndns]
			    end
		            Threads.lock(l)
			    A[ijab,(k-1)*ndns*ndns*nvir + ija]=(Ckib*F_MO[ib+ndns,k+ndns])
		            Threads.unlock(l)
		            Threads.lock(l)
			    A[(k-1)*ndns*ndns*nvir + ija,ijab]=(Ckib*F_MO[ib+ndns,k+ndns])
		            Threads.unlock(l)
			end
		    end

                end
            end
        end
    end
    return A
end

function ext_koopmans(p,elag,n)

    elag_small = view(elag,1:p.nbf5,1:p.nbf5)
    n_sqrt = 1 ./ sqrt.(n)
    @tullio nu[q,p] := -elag_small[q,p]*n_sqrt[q]*n_sqrt[p]

    println("")
    println("---------------------------")
    println(" Extended Koopmans Theorem ")
    println("   Ionization Potentials   ")
    println("---------------------------")

    eigval, eigvec = eigen(nu)
    println(" OM        (eV)")
    println("---------------------------")
    for (i,val) in enumerate(reverse(eigval))
        @printf(" %3i       %7.3f\n",i,real(val)*27.2114)
    end

    println("")
    @printf("EKT IP: %7.3f eV",real(eigval[1])*27.2114)
    println("")
end

function mulliken_pop(bset,p,n,C,S)

    C_nbf5 = view(C,1:p.nbf,1:p.nbf5)

    @tullio nPS[m] := 2*n[i]*C_nbf5[m,i]*C_nbf5[l,i]*S[l,m]

    pop = zeros(p.natoms)

    ifun = 1
    for iatom in 1:p.natoms
        for ibasis in 1:size(bset.basis[iatom])[1]
            nfun = 2*bset.basis[iatom][ibasis].l + 1
            pop[iatom] += sum(nPS[ifun:ifun+nfun-1])
            ifun += nfun
        end
    end

    println("")
    println("---------------------------------")
    println("  Mulliken Population Analysis   ")
    println("---------------------------------")
    println(" Idx  Atom   Population   Charge ")
    println("---------------------------------")
    for iatom in 1:p.natoms
	symbol = Z_to_symbol(bset.atoms[iatom].Z)
        @printf("%3i    %2s    %5.2f      %5.2f\n",iatom, symbol, pop[iatom], bset.atoms[iatom].Z-pop[iatom])
    end

end

function lowdin_pop(bset,p,n,C,S)

    evals,evecs = eigen(S)
    for i in 1:size(evals)[1]
        if (evals[i]<0.0)
	    evals[i] = 0.0
        else
	    evals[i] = sqrt(evals[i])
        end
    end
    S_12 = evecs*Diagonal(evals)*evecs'

    C_nbf5 = view(C,1:p.nbf,1:p.nbf5)
    @tullio tmp[s,l] := 2*S_12[s,m]*n[i]*C_nbf5[m,i]*C_nbf5[l,i]#*Sp_12[l,s]
    @tullio S_12nPS_12[s] := tmp[s,l]*S_12[l,s]

    pop = zeros(p.natoms)

    ifun = 1
    for iatom in 1:p.natoms
	for ibasis in 1:size(bset.basis[iatom])[1]
	    nfun = 2*bset.basis[iatom][ibasis].l + 1
	    pop[iatom] += sum(S_12nPS_12[ifun:ifun+nfun-1])
	    ifun += nfun
        end
    end

    println("")
    println("---------------------------------")
    println("   Lowdin Population Analysis    ")
    println("---------------------------------")
    println(" Idx  Atom   Population   Charge ")
    println("---------------------------------")
    for iatom in 1:p.natoms
	symbol = Z_to_symbol(bset.atoms[iatom].Z)
        @printf("%3i    %2s    %5.2f      %5.2f\n",iatom, symbol, pop[iatom], bset.atoms[iatom].Z-pop[iatom])
    end
end

function M_diagnostic(p,n)

    m_vals = 2*n
    if(p.HighSpin)
        m_vals[p.nbeta:p.nalpha] = n[p.nbeta:p.nalpha]
    end

    m_diagnostic = 0

    m_vals[p.no1+1:p.nbeta] = 2.0 .- m_vals[p.no1+1:p.nbeta]
    m_diagnostic += maximum(m_vals[p.no1+1:p.nbeta])

    #if(p.nsoc!=0): #This is always zero
    #    m_vals[p.nbeta+1:p.nalpha] = 1.0 - m_vals[p.nbeta+1:p.nalpha]
    #    m_diagnostic += max(m_vals[p.nbeta+1:p.nalpha])

    m_vals[p.nalpha+1:p.nbf5] = m_vals[p.nalpha+1:p.nbf5] .- 0.0
    m_diagnostic += maximum(m_vals[p.nalpha+1:p.nbf5])

    m_diagnostic = 0.5 * m_diagnostic

    println("")
    println("---------------------------------")
    @printf("   M Diagnostic: %4.2f\n",m_diagnostic)
    println("---------------------------------")
    println("")

end

function mbpt(n,C,H,I,b_mnl,E_nuc,E_elec,p)

    t1 = time()

    println(" MBPT")
    println("=========")

    nab = p.nvir*p.nalpha
    last_coup = p.nalpha + p.ncwo*(p.ndoc)

    @printf("Number of orbitals        (NBASIS) = %d\n",p.nbf)
    @printf("Number of frozen pairs       (NFR) = %d\n",p.no1)
    @printf("Number of occupied orbs.     (NOC) = %d\n",p.nalpha)
    @printf("Number of virtual orbs.     (NVIR) = %d\n",p.nvir)
    @printf("Size of A+B and A-B (NAB=NOCxNVIR) = %d\n",nab)
    println("")

    occ = zeros(p.nbf)
    occ[1:p.nbf5] = n

    println(" ....Building F_MO")
    EHFL, F_MO = build_F_MO(C,H,I,b_mnl,p)

    println(" ....Transforming ERIs mnsl->pqrt")

    if p.RI
        @tullio tmp[p,n,l] := C[m,p] * b_mnl[m,n,l]
        @tullio pql[p,q,l] := C[n,q] * tmp[p,n,l]
    else
        @tullio tmp[p,n,s,l] := C[m,p] * I[m,n,s,l]
        @tullio tmp2[p,q,s,l] := C[n,q] * tmp[p,n,s,l]
        @tullio tmp3[p,q,r,l] := C[s,r] * tmp2[p,q,s,l]
        @tullio pqrt[p,q,r,t] := C[l,t] * tmp3[p,q,r,l]
    end

    # Eq. (17) and (18)
    Cintra = 1 .- (1 .- abs.(1 .- 2 * occ)).^2
    Cinter = abs.(1 .-2 * occ).^2
    Cinter[1:p.nalpha] .= 1.0

    # Eq. (19) Canonicalization of NOs Step 1
    println(" ....Attenuating F_MO")
    F_MO_at = F_MO_attenuated(F_MO,Cintra,Cinter,p.no1,p.nalpha,p.ndoc,p.nsoc,p.ndns,p.ncwo,p.nbf5,p.nbf)
    
    # Eq. (20)
    println(" ....Attenuating pqrt")
    if p.RI
        pql_at_intra, pql_at_inter = ERIS_RI_attenuated(pql,Cintra,Cinter,p.no1,p.ndoc,p.nsoc,p.ndns,p.ncwo,p.nbf5,p.nbf,p.nbfaux)
    else
        pqrt_at = ERIS_attenuated(pqrt,Cintra,Cinter,p.no1,p.ndoc,p.nsoc,p.ndns,p.ncwo,p.nbf5,p.nbf)
    end

    # Canonicalization of NOs Step 2 and 3
    eig,C_can = eigen(Symmetric(F_MO_at))

    if(p.RI)
        @tullio tmp[p,n,l] := C_can[m,p] * pql[m,n,l]
        @tullio pql[p,q,l] := C_can[n,q] * tmp[p,n,l]

        @tullio tmp[p,n,l] := C_can[m,p] * pql_at_intra[m,n,l]
        @tullio pql_at_intra[p,q,l] := C_can[n,q] * tmp[p,n,l]
        @tullio tmp[p,n,l] := C_can[m,p] * pql_at_inter[m,n,l]
        @tullio pql_at_inter[p,q,l] := C_can[n,q] * tmp[p,n,l]
    else
        @tullio tmp[p,n,s,l] := C_can[m,p] * pqrt[m,n,s,l]
        @tullio tmp2[p,q,s,l] := C_can[n,q] * tmp[p,n,s,l]
        @tullio tmp3[p,q,r,l] := C_can[s,r] * tmp2[p,q,s,l]
        @tullio pqrt[p,q,r,t] := C_can[l,t] * tmp3[p,q,r,l]

        @tullio tmp[p,n,s,l] := C_can[m,p] * pqrt_at[m,n,s,l]
        @tullio tmp2[p,q,s,l] := C_can[n,q] * tmp[p,n,s,l]
        @tullio tmp3[p,q,r,l] := C_can[s,r] * tmp2[p,q,s,l]
        @tullio pqrt_at[p,q,r,t] := C_can[l,t] * tmp3[p,q,r,l]
    end

    println("List of qp-orbital energies (a.u.) and occ numbers used")
    println()
    mu = (eig[p.nalpha] + eig[p.nalpha+1])/2
    eig = eig .- mu
    for i in 1:p.nbf
        @printf(" %8.6f %5.3f\n",eig[i],occ[i])
    end
    @printf("Chemical potential used for qp-orbital energies: %6.4f\n",mu)
    println("")

    if p.RI
        EcMP2 = mp2_RI_eq(eig,pql,pql_at_intra,pql_at_inter,p.no1,p.ndoc,p.nsoc,p.ndns,p.ncwo,p.nbf5,p.nbeta,p.nalpha,p.nbf)
    else
        EcMP2 = mp2_eq(eig,pqrt,pqrt_at,p.nbeta,p.nalpha,p.nbf)
    end

    ECndHF,ECndl = ECorrNonDyn(n,C,H,I,b_mnl,p)

    ESD = EHFL+E_nuc+ECndHF
    ESDc = ESD + ECndl
    EPNOF = E_elec + E_nuc

    @printf(" E(SD)                = %f\n",ESD)
    @printf(" E(SD+ND)             = %f\n",ESDc)
    @printf(" E(PNOFi)             = %f\n",EPNOF)
    @printf("\n")
    @printf(" Ec(ND)               = %f\n",ECndl)
#    print(" Ec(RPA-FURCHE)       = {: f}".format(EcRPA))
#    print(" Ec(RPA)              = {: f}".format(iEcRPA))
#    print(" Ec(AC-SOSEX)         = {: f}".format(iEcSOSEX))
#    print(" Ec(RPA+AC-SOSEX)     = {: f}".format(iEcRPASOS))
#    print(" Ec(GW@GM)            = {: f}".format	(EcGoWo))
#    print(" Ec(SOSEX@GM)         = {: f}".format(EcGMSOS))
#    print(" Ec(GW@GM+SOSEX@GM)   = {: f}".format(EcGoWoSOS))
    @printf(" Ec(MP2)              = %f\n",EcMP2)
#    print(" Ec(CCSD)             = {: f}".format(EcCCSD))
#    print("")
#    print(" E(RPA-FURCHE)       = {: f}".format(ESD + EcRPA))
#    print(" E(RPA)              = {: f}".format(ESD + iEcRPA))
#    print(" E(RPA+AC-SOSEX)     = {: f}".format(ESD + iEcRPASOS))
#    print(" E(GW@GM)            = {: f}".format(ESD + EcGoWo))
#    print(" E(SOSEX@GM)         = {: f}".format(ESD + EcGMSOS))
#    print(" E(GW@GM+SOSEX@GM)   = {: f}".format(ESD + EcGoWoSOS))
    @printf(" E(MP2)               = %f\n",ESD + EcMP2)
#    print(" E(CCSD)             = {: f}".format(ESD + EcCCSD))
#    print("")
#    print(" E(NOF-c-RPA-FURCHE)       = {: f}".format(ESDc + EcRPA))
#    print(" E(NOF-c-RPA)              = {: f}".format(ESDc + iEcRPA))
#    print(" E(NOF-c-RPA+AC+SOSEX)     = {: f}".format(ESDc + iEcRPASOS))
#    print(" E(NOF-c-GW@GM)            = {: f}".format(ESDc + EcGoWo))
#    print(" E(NOF-c-SOSEX@GM)         = {: f}".format(ESDc + EcGMSOS))
#    print(" E(NOF-c-GW@GM+SOSEX@GM)   = {: f}".format(ESDc + EcGoWoSOS))
    @printf(" E(NOF-c-MP2)         = %f\n",ESDc + EcMP2)
#    print(" E(NOF-c-CCSD)             = {: f}".format(ESDc + EcCCSD))


end

function build_F_MO(C,H,I,b_mnl,p)

    C_beta = C[:,1:p.nbeta]
    D = C_beta * C_beta'
    if(p.MSpin==0)
        if(p.nsoc>0)
            C_alpha = C[:,p.nbeta+1:p.nalpha]
	    Dalpha = C_alpha * C_alpha'
            D = D + 0.5*Dalpha
        end

	if(p.gpu)
	    D = CuArray(D)
            @tensor J[m,n] := D[l,s] * I[m,n,s,l]
            @tensor K[m,s] := D[n,l] * I[m,n,s,l]
	    J,K,D = Array(J), Array(K), Array(D)
        else
	    if p.RI
		@tullio X[l] := b_mnl[m,n,l] * D[m,n]
                @tullio J[m,n] := b_mnl[m,n,l] * X[l]
		@tullio tmp[m,s,l] := b_mnl[m,n,l] * D[s,n]
		@tullio K[m,s] := b_mnl[m,n,l] * tmp[s,n,l]
	    else
                @tullio J[m,n] := D[l,s] * I[m,n,s,l]
                @tullio K[m,s] := D[n,l] * I[m,n,s,l]
            end
	end
        F = H + 2*J - K
        EHFL = tr( D * H + D * F )
    end

    F_MO = C' * F * C
    eig = diag(F_MO)

    return EHFL, F_MO

end

function F_MO_attenuated(F_MO,Cintra,Cinter,no1,nalpha,ndoc,nsoc,ndns,ncwo,nbf5,nbf)

    F_MO_at = zeros(nbf,nbf)

    subspaces = zeros(nbf)
    for i in 1:no1
        subspaces[i] = i
    end
    for i in 1:ndoc
        subspaces[no1+i] = no1+i
        ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        subspaces[ll:ul] .= no1+i
    end
    for i in 1:nsoc
        subspaces[no1+ndoc+i] = no1+ndoc+i
    end
    subspaces[nbf5+1:end] .= -1

    for p in 1:nbf
        for q in 1:nbf
            if(p != q)
                if(subspaces[p]==subspaces[q])
                    F_MO_at[p,q] = F_MO[p,q]*Cintra[p]*Cintra[q]
                else
                    F_MO_at[p,q] = F_MO[p,q]*Cinter[p]*Cinter[q]
		end
            else
                F_MO_at[p,q] = F_MO[p,q]
	    end
	end
    end
    F_MO_at[1:nalpha,nalpha+1:nbf] .= 0.0
    F_MO_at[nalpha+1:nbf,1:nalpha] .= 0.0

    return F_MO_at

end

function ERIS_attenuated(pqrt,Cintra,Cinter,no1,ndoc,nsoc,ndns,ncwo,nbf5,nbf)

    subspaces = zeros(nbf)
    for i in 1:no1
        subspaces[i] = i
    end
    for i in 1:ndoc
        subspaces[no1+i] = no1+i
        ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        subspaces[ll:ul] .= no1+i
    end
    for i in 1:nsoc
        subspaces[no1+ndoc+i] = no1+ndoc+i
    end
    subspaces[nbf5+1:end] .= -1

    pqrt_at = zeros(nbf,nbf,nbf,nbf)
    for p in 1:nbf
        for q in 1:nbf
            for r in 1:nbf
                for t in 1:nbf
                    if(subspaces[p]==subspaces[q] && subspaces[p]==subspaces[r] && subspaces[p]==subspaces[t] && subspaces[p]!=-1)
                        pqrt_at[p,q,r,t] = pqrt[p,q,r,t] * Cintra[p]*Cintra[q]*Cintra[r]*Cintra[t]
                    else
                        pqrt_at[p,q,r,t] = pqrt[p,q,r,t] * Cinter[p]*Cinter[q]*Cinter[r]*Cinter[t]
                    end
                end
            end
        end
    end
    return pqrt_at
end

function ERIS_RI_attenuated(pql,Cintra,Cinter,no1,ndoc,nsoc,ndns,ncwo,nbf5,nbf,nbfaux)

    subspaces = zeros(nbf)
    for i in 1:no1
        subspaces[i] = i
    end
    for i in 1:ndoc
        subspaces[no1+i] = no1+i
        ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        subspaces[ll:ul] .= no1+i
    end
    for i in 1:nsoc
        subspaces[no1+ndoc+i] = no1+ndoc+i
    end
    subspaces[nbf5+1:end] .= -1

    pql_at_intra = zeros(nbf,nbf,nbfaux)
    pql_at_inter = zeros(nbf,nbf,nbfaux)
    for p in 1:nbf
        for q in 1:nbf
            for l in 1:nbfaux
                pql_at_intra[p,q,l] = pql[p,q,l] * Cintra[p]*Cintra[q]
                pql_at_inter[p,q,l] = pql[p,q,l] * Cinter[p]*Cinter[q]
            end
        end
    end
    return pql_at_intra, pql_at_inter
end

function mp2_eq(eig,pqrt,pqrt_at,nbeta,nalpha,nbf)

    EcMP2 = 0
    for a in nalpha+1:nbf
        for b in nalpha+1:nbf
            for i in 1:nbeta
                for j in 1:nbeta
                    EcMP2 += pqrt[i,a,j,b]*(2*pqrt_at[i,a,j,b]-pqrt_at[i,b,j,a])/(eig[i]+eig[j]-eig[a]-eig[b]+1e-10)
                end
                for j in nbeta+1:nalpha
                    EcMP2 += pqrt[i,a,j,b]*(pqrt_at[i,a,j,b]-0.5*pqrt_at[i,b,j,a])/(eig[i]+eig[j]-eig[a]-eig[b]+1e-10)
                end
            end
            for i in nbeta+1:nalpha
                for j in 1:nbeta
                    EcMP2 += pqrt[i,a,j,b]*(pqrt_at[i,a,j,b]-0.5*pqrt_at[i,b,j,a])/(eig[i]+eig[j]-eig[a]-eig[b]+1e-10)
                end
                for j in nbeta+1:nalpha
                    EcMP2 += 0.5*pqrt[i,a,j,b]*(pqrt_at[i,a,j,b]-0.5*pqrt_at[i,b,j,a])/(eig[i]+eig[j]-eig[a]-eig[b]+1e-10)
                end
	    end
        end
    end

    return EcMP2
end

function mp2_RI_eq(eig,pql,pql_at_intra,pql_at_inter,no1,ndoc,nsoc,ndns,ncwo,nbf5,nbeta,nalpha,nbf)


    subspaces = zeros(nbf)
    for i in 1:no1
        subspaces[i] = i
    end
    for i in 1:ndoc
        subspaces[no1+i] = no1+i
        ll = no1 + ndns + ncwo*(ndoc - i) + 1
        ul = no1 + ndns + ncwo*(ndoc - i + 1)
        subspaces[ll:ul] .= no1+i
    end
    for i in 1:nsoc
        subspaces[no1+ndoc+i] = no1+ndoc+i
    end
    subspaces[nbf5+1:end] .= -1

    EcMP2 = 0
    for i in 1:nbeta
	ial = pql[i,nalpha+1:nbf,1:end]
	ial_at_intra = pql_at_intra[i,1:nbf,1:end]
	ial_at_inter = pql_at_inter[i,1:nbf,1:end]
	ei = eig[i]
        for j in 1:nbeta
	    jbl = pql[j,nalpha+1:nbf,1:end]
	    jbl_at_intra = pql_at_intra[j,1:nbf,1:end]
	    jbl_at_inter = pql_at_inter[j,1:nbf,1:end]
	    ej = eig[j]

	    @tullio I[a,b] := ial[a,l] * jbl[b,l]
	    e = eig[nalpha+1:nbf]
            I_at = zeros(nbf-nalpha,nbf-nalpha)
	    for a in nalpha+1:nbf
	        for b in nalpha+1:nbf
		    if(subspaces[i]==subspaces[a] && subspaces[i]==subspaces[b] && subspaces[i]==subspaces[j] && subspaces[i]!=-1)
			ial_tmp = ial_at_intra[a,1:end]
			jbl_tmp = jbl_at_intra[b,1:end]
		    else
			ial_tmp = ial_at_inter[a,1:end]
			jbl_tmp = jbl_at_inter[b,1:end]
	            end
		    @tullio tmp = ial_tmp[l] * jbl_tmp[l]
		    I_at[a-nalpha,b-nalpha] = tmp
	        end
	    end
	    @tullio EcMP2 += I[a,b]*(2*I_at[a,b] - I_at[b,a])/(ei+ej-e[a]-e[b])
#            EcMP2 += pqrt[i,a,j,b]*(2*pqrt_at[i,a,j,b]-pqrt_at[i,b,j,a])/(eig[i]+eig[j]-eig[a]-eig[b]+1e-10)
        end
        for j in nbeta+1:nalpha
            jbl = pql[j,nalpha+1:nbf,1:end]
            jbl_at_intra = pql_at_intra[j,1:nbf,1:end]
            jbl_at_inter = pql_at_inter[j,1:nbf,1:end]
            ej = eig[j]

            @tullio I[a,b] := ial[a,l] * jbl[b,l]
            e = eig[nalpha+1:nbf]
            I_at = zeros(nbf-nalpha,nbf-nalpha)
            for a in nalpha+1:nbf
                for b in nalpha+1:nbf
                    if(subspaces[i]==subspaces[a] && subspaces[i]==subspaces[b] && subspaces[i]==subspaces[j] && subspaces[i]!=-1)
                        ial_tmp = ial_at_intra[a,1:end]
                        jbl_tmp = jbl_at_intra[b,1:end]
                    else
                        ial_tmp = ial_at_inter[a,1:end]
                        jbl_tmp = jbl_at_inter[b,1:end]
                    end
		    @tullio tmp = ial_tmp[l] * jbl_tmp[l]
		    I_at[a-nalpha,b-nalpha] = tmp
                end
            end
	    @tullio EcMP2 += I[a,b]*(I_at[a,b] - 0.5*I_at[b,a])/(ei+ej-e[a]-e[b])
#            EcMP2 += pqrt[i,a,j,b]*(pqrt_at[i,a,j,b]-0.5*pqrt_at[i,b,j,a])/(eig[i]+eig[j]-eig[a]-eig[b]+1e-10)
        end
    end
    for i in nbeta+1:nalpha
	ial = pql[i,nalpha+1:nbf,1:end]
	ial_at_intra = pql_at_intra[i,1:nbf,1:end]
	ial_at_inter = pql_at_inter[i,1:nbf,1:end]
	ei = eig[i]
        for j in 1:nbeta
            jbl = pql[j,nalpha+1:nbf,1:end]
            jbl_at_intra = pql_at_intra[j,1:nbf,1:end]
            jbl_at_inter = pql_at_inter[j,1:nbf,1:end]
            ej = eig[j]

            @tullio I[a,b] := ial[a,l] * jbl[b,l]
            e = eig[nalpha+1:nbf]
            I_at = zeros(nbf-nalpha,nbf-nalpha)
            for a in nalpha+1:nbf
                for b in nalpha+1:nbf
                    if(subspaces[i]==subspaces[a] && subspaces[i]==subspaces[b] && subspaces[i]==subspaces[j] && subspaces[i]!=-1)
                        ial_tmp = ial_at_intra[a,1:end]
                        jbl_tmp = jbl_at_intra[b,1:end]
                    else
                        ial_tmp = ial_at_inter[a,1:end]
                        jbl_tmp = jbl_at_inter[b,1:end]
                    end
		    @tullio tmp = ial_tmp[l] * jbl_tmp[l]
		    I_at[a-nalpha,b-nalpha] = tmp
                end
            end
	    @tullio EcMP2 += I[a,b]*(I_at[a,b] - 0.5*I_at[b,a])/(ei+ej-e[a]-e[b])
#            EcMP2 += pqrt[i,a,j,b]*(pqrt_at[i,a,j,b]-0.5*pqrt_at[i,b,j,a])/(eig[i]+eig[j]-eig[a]-eig[b]+1e-10)
        end
        for j in nbeta+1:nalpha
            jbl = pql[j,nalpha+1:nbf,1:end]
            jbl_at_intra = pql_at_intra[j,1:nbf,1:end]
            jbl_at_inter = pql_at_inter[j,1:nbf,1:end]
            ej = eig[j]

            @tullio I[a,b] := ial[a,l] * jbl[b,l]
            e = eig[nalpha+1:nbf]
            I_at = zeros(nbf-nalpha,nbf-nalpha)
            for a in nalpha+1:nbf
                for b in nalpha+1:nbf
                    if(subspaces[i]==subspaces[a] && subspaces[i]==subspaces[b] && subspaces[i]==subspaces[j] && subspaces[i]!=-1)
                        ial_tmp = ial_at_intra[a,1:end]
                        jbl_tmp = jbl_at_intra[b,1:end]
                    else
                        ial_tmp = ial_at_inter[a,1:end]
                        jbl_tmp = jbl_at_inter[b,1:end]
                    end
		    @tullio tmp = ial_tmp[l] * jbl_tmp[l]
		    I_at[a-nalpha,b-nalpha] = tmp
                end
            end
            @tullio EcMP2 += 0.5*I[a,b]*(I_at[a,b] - 0.5*I_at[b,a])/(ei+ej-e[a]-e[b])		
#            EcMP2 += 0.5*pqrt[i,a,j,b]*(pqrt_at[i,a,j,b]-0.5*pqrt_at[i,b,j,a])/(eig[i]+eig[j]-eig[a]-eig[b]+1e-10)
        end
    end

    return EcMP2
end

function ECorrNonDyn(n,C,H,I,b_mnl,p)

    fi = 2 .* n .* (1 .- n)
    CK12nd = fi .* fi'

    beta = sqrt.((1 .- abs.(1 .-2 .* n)) .*n )

    for l in 1:p.ndoc
        ll = p.no1 + p.ndns + p.ncwo*(p.ndoc - l) + 1
        ul = p.no1 + p.ndns + p.ncwo*(p.ndoc - l + 1)
        CK12nd[p.no1+l,ll:ul] .= beta[p.no1+l] .* beta[ll:ul]
        CK12nd[ll:ul,p.no1+l] .= beta[ll:ul] .* beta[p.no1+l]
        CK12nd[ll:ul,ll:ul] .= -beta[ll:ul] .* beta[ll:ul]'
    end

    #C^K KMO
    J_MO,K_MO,H_core = computeJKH_MO(C,H,I,b_mnl,p)

    ECndHF = 0
    ECndl = 0
    if (p.MSpin==0)
	
        ECndHF = - dot(diag(CK12nd[p.nbeta+1:p.nalpha,p.nbeta+1:p.nalpha]),diag(K_MO[p.nbeta+1:p.nalpha,p.nbeta+1:p.nalpha]))
	CK12nd[diagind(CK12nd)] .= 0
	ECndl -= sum(CK12nd .* K_MO')
    end

    return ECndHF,ECndl
end
