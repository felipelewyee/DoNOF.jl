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
        @printf(" %3i       %7.3f\n",i,val*27.2114)
    end

    println("")
    @printf("EKT IP: %7.3f eV",eigval[1]*27.2114)
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

    m_vals[p.no1:p.nbeta] = 2.0 .- m_vals[p.no1:p.nbeta]
    m_diagnostic += maximum(m_vals[p.no1:p.nbeta])

    #if(p.nsoc!=0): #This is always zero
    #    m_vals[p.nbeta:p.nalpha] = 1.0 - m_vals[p.nbeta:p.nalpha]
    #    m_diagnostic += max(m_vals[p.nbeta:p.nalpha]) 

    m_vals[p.nalpha:p.nbf5] = m_vals[p.nalpha:p.nbf5] .- 0.0
    m_diagnostic += maximum(m_vals[p.nalpha:p.nbf5])

    println("")
    println("---------------------------------")
    @printf("   M Diagnostic: %4.2f\n",m_diagnostic)
    println("---------------------------------")
    println("")

end
