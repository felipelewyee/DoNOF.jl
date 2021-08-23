function nofmp2(n,C,H,I,b_mnl,E_nuc,p)

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
        J,K = JK_HF_Full(D,I,p)
        #J,K = computeJK_HF(D,I,b_mnl,p)
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

    #println(EHFL)

    @tullio F_MO[i,j] := vec[k,i]*F[k,l]*vec[l,j]
    #F_MO = np.matmul(np.matmul(np.transpose(vec),F),vec)

    F_MO_act = view(F_MO,1:p.nbf-p.no1,1:p.nbf-p.no1)
    @tullio eig[i] := F_MO_act[i,i]
    #eig = np.einsum("ii->i",F_MO[:p.nbf-p.no1,:p.nbf-p.no1])

    iajb = iajb_Full_jit(C,I,p.no1,p.nalpha,p.nbf,p.nbf5)
    #iajb = compute_iajb(C,I,b_mnl,p)
    
    FI1 = ones(p.nbf-p.no1)
    FI2 = ones(p.nbf-p.no1)

    FI1[1:p.nbf5-p.no1] = 1.0 .- (1.0 .- abs.(1.0 .-2*occ[1:p.nbf5-p.no1])).^2

    FI2[p.nalpha-p.no1:p.nbf5-p.no1] = abs.(1.0 .-2*occ[p.nalpha-p.no1:p.nbf5-p.no1]).^2

    Tijab = CalTijab(iajb,F_MO,eig,FI1,FI2,p)
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
            CK12nd_alpha_beta = view(CK12nd_alpha_beta,p.nbeta:p.nalpha,p.nbeta:p.nalpha)
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

function CalTijab(iajb,F_MO,eig,FI1,FI2,p)

    println("Starting CalTijab")

    B = build_B(iajb,FI1,FI2,p.ndoc,p.ndns,p.nvir,p.ncwo)
    println("....B vector Computed")

    Tijab = Tijab_guess(iajb,eig,p.ndoc,p.ndns,p.nvir)
    println("....Tijab Guess Computed")

    #A_CSR = csr_matrix(build_A(F_MO,FI1,FI2,p.no1,p.ndoc,p.ndns,p.nvir,p.ncwo,p.nbf))
    A = build_A(F_MO,FI1,FI2,p.no1,p.ndoc,p.ndns,p.nvir,p.ncwo,p.nbf)
    println("A matrix built")
    #"has {}/{} elements with Tol = {}".format(len(A),p.nvir**4*p.ndoc**4,1e-10)    
    #Tijab = solve_Tijab(A_CSR,B,Tijab,p)

    Tijab = B/A
    println("Tijab found")
    #res = root(build_R, Tijab, args=(B,F_MO,FI1,FI2,p.no1,p.ndoc,p.ndns,p.nvir,p.ncwo,p.nbf),method="krylov")
    #if(res.success):
    #    print("....Tijab found as a Root of R = B - A*Tijab in {} iterations".format(res.nit))
    #else:
    #    print("....WARNING! Tijab NOT FOUND as a Root of R = B - A*Tijab in {} iterations".format(res.nit))
    #    print(res)
    #Tijab = res.x
    #print("")

    return Tijab

end

function build_B(iajb,FI1,FI2,ndoc,ndns,nvir,ncwo)
    B = zeros(ndns^2*nvir^2)
    for i in 1:ndns
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
    return transpose(B)

end

function Tijab_guess(iajb,eig,ndoc,ndns,nvir)
    Tijab = zeros(nvir^2*ndns^2)
    for ia in 1:nvir
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
        npair[ll:ul] = i
    end

    Bp = zeros(ndns^2*nvir^2)

    for ib in 1:nvir
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
		        end
			Bp[ijab] += (- Cki*F_MO[i,k])*T[(k-1)+jab]
		    end
                    
                    for k in i+1:ndns
			if abs(F_MO[i,k])>1e-10
                            Cki = FI2[k]*FI2[i]
		        end
			Bp[ijab] += (- Cki*F_MO[i,k])*T[(k-1)+jab]
		    end

                    for k in 1:j-1
                        if abs(F_MO[j,k])>1e-10
                            Ckj = FI2[k]*FI2[j]
		        end
			Bp[ijab] += (- Ckj*F_MO[j,k])*T[(k-1)*ndns+iab]
		    end
                    for k in j+1:ndns
                        if abs(F_MO[j,k])>1e-10
                            Ckj = FI2[k]*FI2[j]
		        end
    		        Bp[ijab] += (- Ckj*F_MO[j,k])*T[(k-1)*ndns+iab]
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
    return R
end



function build_A(F_MO,FI1,FI2,no1,ndoc,ndns,nvir,ncwo,nbf)
    npair = zeros(nvir)
    for i in 1:ndoc
        ll = ncwo*(ndoc - i) + 1
        ul = ncwo*(ndoc - i + 1)
        npair[ll:ul] .= i
    end
    A = spzeros(ndns^2*nvir^2,ndns^2*nvir^2)
    #IROW = np.empty((2*ndns**2*nvir**2*(nbf-no1)),dtype=np.int32)
    #ICOL = np.empty((2*ndns**2*nvir**2*(nbf-no1)),dtype=np.int32)

    #nnz = -1
    for ib in 1:nvir
        for ia in 1:nvir
            for j in 1:ndns
                for i in 1:ndns
                    #print(nnz)
                    jab =     (j-1)*ndns + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir
                    iab = i              + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir
                    ijb = i + (j-1)*ndns                    + (ib-1)*ndns*ndns*nvir
                    ija = i + (j-1)*ndns + (ia-1)*ndns*ndns
                    ijab= i + (j-1)*ndns + (ia-1)*ndns*ndns + (ib-1)*ndns*ndns*nvir

                    #nnz = nnz + 1
                    A[ijab,i + jab] = (F_MO[ia+ndns,ia+ndns] + F_MO[ib+ndns,ib+ndns] - F_MO[i,i] - F_MO[j,j])
                    #IROW[nnz] = (ijab)
                    #ICOL[nnz] = (i + jab)

                    for k in 1:i-1
                        if abs(F_MO[i,k])>1e-10
                            Cki = FI2[k]*FI2[i]
                            #nnz += 1
                            A[ijab,k + jab]=(- Cki*F_MO[i,k])
                            #IROW[nnz]=(ijab)
                            #ICOL[nnz]=(k + jab)
                            #nnz += 1
                            A[k+jab,ijab]=(- Cki*F_MO[i,k])
                            #ICOL[nnz]=(ijab)
                            #IROW[nnz]=(k + jab)
			end
		    end
                    for k in 1:j-1
                        if abs(F_MO[j,k])>1e-10
                            Ckj = FI2[k]*FI2[j]
                            #nnz += 1
			    A[ijab,(k-1)*ndns+iab]=(- Ckj*F_MO[j,k])
                            #IROW[nnz]=(ijab)
                            #ICOL[nnz]=(k*ndns + iab)
                            #nnz += 1
                            #A[nnz]=(- Ckj*F_MO[j,k])
			    A[(k-1)*ndns+iab,ijab]=(- Ckj*F_MO[j,k])
                            #ICOL[nnz]=(ijab)
                            #IROW[nnz]=(k*ndns + iab)
			end
		    end

                    for k in 1:ia-1
                        if abs(F_MO[ia+ndns,k+ndns])>1e-10
                            if npair[k]==npair[ia]
                                Ckia = FI1[k+ndns]*FI1[ia+ndns]
                            else
                                Ckia = FI2[k+ndns]*FI2[ia+ndns]
			    end
                            #nnz += 1
			    A[ijab,(k-1)*ndns*ndns + ijb]=(Ckia*F_MO[ia+ndns,k+ndns])
                            #IROW[nnz]=(ijab)
                            #ICOL[nnz]=(k*ndns*ndns + ijb)
                            #nnz += 1
			    A[(k-1)*ndns*ndns + ijb,ijab]=(Ckia*F_MO[ia+ndns,k+ndns])
                            #ICOL[nnz]=(ijab)
                            #IROW[nnz]=(k*ndns*ndns + ijb)
			end
		    end

                    for k in 1:ib-1
                        if abs(F_MO[ib+ndns,k+ndns])>1e-10
                            if npair[k]==npair[ib]
                                Ckib = FI1[k+ndns]*FI1[ib+ndns]
                            else
                                Ckib = FI2[k+ndns]*FI2[ib+ndns]
			    end
                            #nnz += 1
			    A[ijab,(k-1)*ndns*ndns*nvir + ija]=(Ckib*F_MO[ib+ndns,k+ndns])
                            #IROW[nnz]=(ijab)
                            #ICOL[nnz]=(k*ndns*ndns*nvir + ija)
                            #nnz += 1
                            #A[nnz]=(Ckib*F_MO[ib+ndns,k+ndns])
			    A[(k-1)*ndns*ndns*nvir + ija,ijab]=(Ckib*F_MO[ib+ndns,k+ndns])
                            #ICOL[nnz]=(ijab)
                            #IROW[nnz]=(k*ndns*ndns*nvir + ija)
			end
		    end

    #A = A[:nnz+1]
    #IROW = IROW[:nnz+1]
    #ICOL = ICOL[:nnz+1]

end
end
end
end
    return A#,(IROW,ICOL)
end


