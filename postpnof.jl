module postpnof

using Tullio
include("integrals.jl")

function nofmp2(n,C,H,I,b_mnl,E_nuc,p)

    println(" NOF-MP2")
    println("=========")

    occ = view(n,p.no1+1:p.nbf5)
    vec = view(C,:,p.no1+1:p.nbf)

    D = integrals.computeD_HF(C,I,b_mnl,p)
    if p.MSpin==0
        if p.nsoc>0
            Dalpha = integrals.computeDalpha_HF(C,I,b_mnl,p)
            D = D + 0.5*Dalpha
        end
        J,K = integrals.JK_HF_Full(D,I,p)
        #J,K = integrals.computeJK_HF(D,I,b_mnl,p)
        F = H + 2*J - K
	@tullio DH[i,j] := D[i,k]*H[k,j]
	@tullio DF[i,j] := D[i,k]*F[k,j]
	DHDF = DH + DF
	@tullio EHFL := DHDF[i,i]
        #EHFL = np.trace(np.matmul(D,H)+np.matmul(D,F))
    elseif p.MSpin!=0
        D = integrals.computeD_HF(C,I,b_mnl,p)
        Dalpha = integrals.computeDalpha_HF(C,I,b_mnl,p)
        J,K = integrals.computeJK_HF(D,I,b_mnl,p)
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
            J,K = integrals.computeJK_HF(0.5*Dalpha,I,b_mnl,p)
            Falpha = J - K
	    @tullio DaFa[i,j] := 0.5*Dalpha[i,k]*Falpha[k,j]
	    @tullio EHFL += 2*DaFa[i,i]
            #EHFL = EHFL + 2*np.trace(np.matmul(0.5*Dalpha,Falpha))
            F = F + Falpha
        end
    end

    println(EHFL)

    #F_MO = np.matmul(np.matmul(np.transpose(vec),F),vec)

    #eig = np.einsum("ii->i",F_MO[:p.nbf-p.no1,:p.nbf-p.no1])

    #iajb = integrals.compute_iajb(C,I,b_mnl,p)
    #FI1 = np.ones(p.nbf-p.no1)
    #FI2 = np.ones(p.nbf-p.no1)

    #FI1[:p.nbf5-p.no1] = 1 - (1 - abs(1-2*occ[:p.nbf5-p.no1]))**2

    #FI2[p.nalpha-p.no1:p.nbf5-p.no1] = abs(1-2*occ[p.nalpha-p.no1:p.nbf5-p.no1])**2

    #Tijab = CalTijab(iajb,F_MO,eig,FI1,FI2,p)
    #ECd = 0
    #for k in range(p.nvir):
    #    for l in range(p.nvir):
    #        for i in range(p.ndoc):
    #            for j in range(p.ndoc):
    #                Xijkl = iajb[j,k,i,l]
    #                ijkl = i+j*p.ndns+k*p.ndns*p.ndns+l*p.ndns*p.ndns*p.nvir
    #                ijlk = i+j*p.ndns+l*p.ndns*p.ndns+k*p.ndns*p.ndns*p.nvir
    #                ECd = ECd + Xijkl*(2*Tijab[ijkl]-Tijab[ijlk])
    #            for j in range(p.ndoc,p.ndns):
    #                Xijkl = iajb[j,k,i,l]
    #                ijkl = i+j*p.ndns+k*p.ndns*p.ndns+l*p.ndns*p.ndns*p.nvir
    #                ijlk = i+j*p.ndns+l*p.ndns*p.ndns+k*p.ndns*p.ndns*p.nvir
    #                ECd = ECd + Xijkl*(Tijab[ijkl]-0.5*Tijab[ijlk])
    #        for i in range(p.ndoc,p.ndns):
    #            for j in range(p.ndoc):
    #                Xijkl = iajb[j,k,i,l]
    #                ijkl = i+j*p.ndns+k*p.ndns*p.ndns+l*p.ndns*p.ndns*p.nvir
    #                ijlk = i+j*p.ndns+l*p.ndns*p.ndns+k*p.ndns*p.ndns*p.nvir
    #                ECd = ECd + Xijkl*(Tijab[ijkl]-0.5*Tijab[ijlk])
    #            for j in range(p.ndoc,p.ndns):
    #                Xijkl = iajb[j,k,i,l]
    #                ijkl = i+j*p.ndns+k*p.ndns*p.ndns+l*p.ndns*p.ndns*p.nvir
    #                ijlk = i+j*p.ndns+l*p.ndns*p.ndns+k*p.ndns*p.ndns*p.nvir
    #                if(j!=i):
    #                    ECd = ECd + Xijkl*(Tijab[ijkl]-0.5*Tijab[ijlk])/2

    #fi = 2*n*(1-n)

    #CK12nd = np.outer(fi,fi)

    #beta = np.sqrt((1-abs(1-2*n))*n)

    #for l in range(p.ndoc):
    #    ll = p.no1 + p.ndns + p.ncwo*(p.ndoc - l - 1)
    #    ul = p.no1 + p.ndns + p.ncwo*(p.ndoc - l)
    #    CK12nd[p.no1+l,ll:ul] = beta[p.no1+l]*beta[ll:ul]
    #    CK12nd[ll:ul,p.no1+l] = beta[ll:ul]*beta[p.no1+l]
    #    CK12nd[ll:ul,ll:ul] = -np.outer(beta[ll:ul],beta[ll:ul])

    #C^K KMO
    #J_MO,K_MO,H_core = integrals.computeJKH_MO(C,H,I,b_mnl,p)

    #ECndHF = 0
    #ECndl = 0
    #if (p.MSpin==0):
    #   ECndHF = - np.einsum('ii,ii->',CK12nd[p.nbeta:p.nalpha,p.nbeta:p.nalpha],K_MO[p.nbeta:p.nalpha,p.nbeta:p.nalpha]) # sum_ij
    #   ECndl -= np.einsum('ij,ji->',CK12nd,K_MO) # sum_ij
    #   ECndl += np.einsum('ii,ii->',CK12nd,K_MO) # Quita i=j
    #elif (not p.MSpin==0):
    #   ECndl -= np.einsum('ij,ji->',CK12nd[p.no1:p.nbeta,p.no1:p.nbeta],K_MO[p.no1:p.nbeta,p.no1:p.nbeta]) # sum_ij
    #   ECndl -= np.einsum('ij,ji->',CK12nd[p.no1:p.nbeta,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,p.no1:p.nbeta]) # sum_ij
    #   ECndl -= np.einsum('ij,ji->',CK12nd[p.nalpha:p.nbf5,p.no1:p.nbeta],K_MO[p.no1:p.nbeta,p.nalpha:p.nbf5]) # sum_ij
    #   ECndl -= np.einsum('ij,ji->',CK12nd[p.nalpha:p.nbf5,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5]) # sum_ij
    #   ECndl += np.einsum('ii,ii->',CK12nd[p.no1:p.nbeta,p.no1:p.nbeta],K_MO[p.no1:p.nbeta,p.no1:p.nbeta]) # Quita i=j
    #   ECndl += np.einsum('ii,ii->',CK12nd[p.nalpha:p.nbf5,p.nalpha:p.nbf5],K_MO[p.nalpha:p.nbf5,p.nalpha:p.nbf5]) # Quita i=j

    #print("      Ehfc      = {:f}".format(EHFL+E_nuc+ECndHF))
    #print("")
    #print("      ECd       = {:f}".format(ECd))
    #print("      ECnd      = {:f}".format(ECndl))
    #print("      Ecorre    = {:f}".format(ECd+ECndl))
    #print("      E(NOFMP2) = {:f}".format(EHFL+ECd+ECndl+E_nuc+ECndHF))
    #print("")



end

end
