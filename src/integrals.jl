function compute_integrals(bset,p)

    # Overlap, Kinetics, Potential
    S = overlap(bset)
    T = kinetic(bset)
    V = nuclear(bset)
    H = T + V
    I = Float64[]
    b_mnl = Float64[]
    if (!p.RI)
        # Integrales de Repulsión Electrónica, ERIs (mu nu | sigma lambda)
	I = ERI_2e4c(bset)
    else
        if p.spherical
	    aux = BasisSet(bset.name*"-jkfit",bset.atoms)
        else
	    aux = BasisSet(bset.name*"-jkfit",bset.atoms,spherical=false,lib=:acsint)
        end
        Iaux = ERI_2e3c(bset,aux)
        G = ERI_2e2c(aux)

        evals,evecs = eigen(G)
        sqrtinv = Float64[]
        for i in 1:size(evals)[1]
            if (evals[i]<0.0)
                append!(sqrtinv,0.0)
            else
                append!(sqrtinv,1/sqrt(evals[i]))
            end
        end
        Gmsqrt = evecs*Diagonal(sqrtinv)*evecs'
        @tullio b_mnl[m,n,l] := Iaux[m,n,k]*Gmsqrt[k,l]

        p.nbfaux = size(b_mnl)[3]
    end

    #if(p.gpu)
    #    if (!p.RI)
    #        I = CuArray(I)
    #    else
    #	    if(p.gpu_bits==64)
    #            b_mnl = CuArray(b_mnl)
    #	    else
    #            b_mnl = cu(b_mnl)
    #	    end
    #    end
    #end

    return S,T,V,H,I,b_mnl

end


######################################### J_mn^(j) K_mn^(j) #########################################

function computeJKj(C,I,b_mnl,p)

    #if(p.gpu)
    #    if(p.RI)
    #        J,K = JKj_RI(C,b_mnl,p.nbf,p.nbf5,p.nbfaux)
    #    else
    #        J,K = JKj_Full(CuArray(C),I,p.nbf,p.nbf5)
    #    end
    #    return J,K
    #else
        if(p.RI)
            J,K = JKj_RI(C,b_mnl,p.nbf,p.nbf5,p.nbfaux)
        else
            J,K = JKj_Full(C,I,p.nbf,p.nbf5)
        end
        return J,K
    #end

end

#########################################

function JKj_Full(C,I,nbf,nbf5)

    Cnbf5 = view(C,:,1:nbf5)
    
    @tullio D[i,m,n] := Cnbf5[m,i]*Cnbf5[n,i]
    @tullio J[i,m,n] := D[i,s,l]*I[m,n,s,l]
    @tullio K[i,m,s] := D[i,n,l]*I[m,n,s,l]

    return J,K

end

function JKj_RI(C,b_mnl,nbf,nbf5,nbfaux)

    Cnbf5 = view(C,:,1:nbf5)
    
    #b_transform
    @tullio b_qnl[q,n,l] := Cnbf5[m,q]*b_mnl[m,n,l]
    @tullio b_qql[q,l] := Cnbf5[n,q]*b_qnl[q,n,l]

    #hstarj
    @tullio J[q,m,n] := b_qql[q,l]*b_mnl[m,n,l]

    #hstark
    @tullio K[q,m,n] := b_qnl[q,m,l]*b_qnl[q,n,l]

    return J,K

end

#function JKj_RI(C,b_mnl::CuArray,nbf,nbf5,nbfaux)
#
#    Cnbf5 = C[1:nbf,1:nbf5]
#    Cnbf5 = CuArray{typeof(b_mnl).parameters[1]}(Cnbf5)
#
#    #b_transform
#    @tensor b_qnl[q,n,l] := Cnbf5[m,q]*b_mnl[m,n,l]
#
#    #hstark
#    #@tullio K[q,m,n] := b_qnl[q,m,l]*b_qnl[q,n,l]
#    K = CUDA.zeros(typeof(b_mnl).parameters[1],nbf5,nbf,nbf)
#    for q in 1:nbf5
#        b_q = b_qnl[q,1:nbf,1:nbfaux]
#	K[q,1:nbf,1:nbf] = b_q * b_q'
#        CUDA.unsafe_free!(b_q)
#    end
#
#    #b_qql = dropdims( sum(Cnbf5' .* b_qnl, dims=2), dims=2)
#    @tullio b_qql[q,l] := Cnbf5[n,q]*b_qnl[q,n,l]
#    CUDA.unsafe_free!(b_qnl)
#    CUDA.unsafe_free!(Cnbf5)
#
#    #hstarj
#    @tensor J[q,m,n] := b_qql[q,l]*b_mnl[m,n,l]
#    CUDA.unsafe_free!(b_qql)
#
#    return J,K
#
#end

######################################### J_pq K_pq #########################################

function computeJKH_MO(C,H,I,b_mnl,p)

    #if(p.gpu)
    #    if(p.RI)
    #        J_MO,K_MO,H_core = JKH_MO_RI(CuArray{typeof(b_mnl).parameters[1]}(C),CuArray{typeof(b_mnl).parameters[1]}(H),b_mnl,p.nbf,p.nbf5,p.nbfaux)
    #    else
    #        J_MO,K_MO,H_core = JKH_MO_Full(CuArray(C),CuArray(H),I,p.nbf,p.nbf5)
    #    end
    #    return Array(J_MO),Array(K_MO),Array(H_core)
    # else
        if(p.RI)
             J_MO,K_MO,H_core = JKH_MO_RI(C,H,b_mnl,p.nbf,p.nbf5,p.nbfaux)
         else
             J_MO,K_MO,H_core = JKH_MO_Full(C,H,I,p.nbf,p.nbf5)
         end
         return J_MO,K_MO,H_core
    #  end

end

#########################################

function JKH_MO_RI(C,H,b_mnl,nbf,nbf5,nbfaux)

    Cnbf5 = view(C,:,1:nbf5)

    #denmatj
    @tullio D[i,m,n] := Cnbf5[m,i]*Cnbf5[n,i]

    #b transform
    @tullio b_pnl[p,n,l] := Cnbf5[m,p] * b_mnl[m,n,l]
    @tullio b_pql[p,q,l] := Cnbf5[n,q] * b_pnl[p,n,l]

    #QJMATm
    @tullio J_MO[p,q] := b_pql[p,p,l]*b_pql[q,q,l]

    #QKMATm
    @tullio K_MO[p,q] := b_pql[p,q,l]*b_pql[p,q,l]

    #QHMATm
    @tullio H_core[i] := D[i,m,n]*H[m,n]

    return J_MO,K_MO,H_core
end

#function JKH_MO_RI(C,H,b_mnl::CuArray,nbf,nbf5,nbfaux)
#
#    Cnbf5 = C[1:nbf,1:nbf5]
#
#    #b transform
#    @tensor b_pnl[p,n,l] := Cnbf5[m,p] * b_mnl[m,n,l]
#    @tensor b_pql[p,q,l] := Cnbf5[n,q] * b_pnl[p,n,l]
#    CUDA.unsafe_free!(b_pnl)
#
#    #QJMATm
#    @tullio tmp[p,l] := b_pql[p,p,l]
#    J_MO = tmp * tmp'
#    CUDA.unsafe_free!(tmp)
#    #@tullio J_MO[p,q] := b_pql[p,p,l]*b_pql[q,q,l]
#
#    #QKMATm
#    @tullio K_MO[p,q] := b_pql[p,q,l]*b_pql[p,q,l]
#    #K_MO = dropdims( sum(b_pql .* b_pql, dims=3), dims=3)
#    CUDA.unsafe_free!(b_pql)
#
#    #QHMATm
#    @tensor tmp[m,i] := H[m,n]*Cnbf5[n,i]
#    @tullio H_core[i] := tmp[m,i]*Cnbf5[m,i]
#    CUDA.unsafe_free!(tmp)
#    #@tullio D[i,m,n] := Cnbf5[m,i]*Cnbf5[n,i]
#    #@tensor H_core[i] := D[i,m,n]*H[m,n]
#
#    return J_MO,K_MO,H_core
#end

function JKH_MO_Full(C,H,I,nbf,nbf5)


    Cnbf5 = view(C,:,1:nbf5)

    #denmatj
    @tullio D[i,m,n] := Cnbf5[m,i]*Cnbf5[n,i]

    #QJMATm
    @tullio J[j,m,n] := D[j,s,l]*I[m,n,s,l]
    @tullio J_MO[i,j] := D[i,m,n]*J[j,m,n]

    #QKMATm
    @tullio K[j,m,s] := D[j,n,l]*I[m,n,s,l]
    @tullio K_MO[i,j] := D[i,m,s]*K[j,m,s]

    #QHMATm
    @tullio H_core[i] := D[i,m,n]*H[m,n]

    return J_MO,K_MO,H_core
end

function computeD_HF(C,I,b_mnl,p)

    Cbeta = view(C,:,1:p.nbeta)

    @tullio D[m,n] := Cbeta[m,j]*Cbeta[n,j]

    return D
end

function computeDalpha_HF(C,I,b_mnl,p)

    Calpha = view(C,:,p.nbeta+1:p.nalpha)
    @tullio D[m,n] := Calpha[m,j]*Calpha[n,j]

    return D
end

function computeJK_HF(D,I,p)

    #if(p.gpu)
    #    J,K = JK_HF_Full(CuArray(D),I,p)
    #    return Array(J),Array(K)
    #else
        J,K = JK_HF_Full(D,I,p)
        return J,K
    #end


end

function JK_HF_Full(D,I,p)

    #denmatj
    @tullio J[m,n] := D[l,s]*I[m,n,s,l]
    @tullio K[m,s] := D[n,l]*I[m,n,s,l]

    return J,K

end

function compute_iajb(C,I,p)

    #if(p.gpu)
    #    iajb = iajb_Full(CuArray(C),I,p.no1,p.nalpha,p.nbf,p.nbf5)
    #	return Array(iajb)
    #else
        iajb = iajb_Full(C,I,p.no1,p.nalpha,p.nbf,p.nbf5)
        return iajb
    #end

end

function iajb_Full(C,I,no1,nalpha,nbf,nbf5)

    Cocc = view(C,:,no1+1:nalpha)
    Ccwo = view(C,:,nalpha+1:nbf)

    @tullio iajb[i,a,j,b] := ((Ccwo[n,a]*((Cocc[m,i]*I[m,n,s,l])*Cocc[s,j]))*Ccwo[l,b])

    return iajb

end
