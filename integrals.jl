module integrals

using Tullio
using LoopVectorization
using PyCall
using CUDA, CUDAKernels, KernelAbstractions # Now defined with a GPU version:

function compute_integrals(wfn,mol,p)

    psi4 = pyimport_conda("psi4", "psi4")
    np = pyimport_conda("numpy","psi4")

    # Integrador
    mints = psi4.core.MintsHelper(wfn.basisset())
    
    # Overlap, Kinetics, Potential
    S = copy(np.asarray(mints.ao_overlap()))
    T = copy(np.asarray(mints.ao_kinetic()))
    V = copy(np.asarray(mints.ao_potential()))
    H = T + V
    I = []
    b_mnl = []
    if (!p.RI)
        # Integrales de Repulsión Electrónica, ERIs (mu nu | sigma lambda)
        I = copy(np.asarray(mints.ao_eri()))
    else
    
        orb = wfn.basisset()
        aux = psi4.core.BasisSet.build(mol, "DF_BASIS_SCF", "", "JKFIT", orb.blend())
        zero_bas = psi4.core.BasisSet.zero_ao_basis_set()
    
        Ppq = mints.ao_eri(orb, orb, aux, zero_bas)
    
        metric = mints.ao_eri(aux, zero_bas, aux, zero_bas)
        metric.power(-0.5, 1.e-14)
        p.nbfaux = metric.shape[1]
    
        Ppq = copy(np.squeeze(Ppq))
        metric = copy(np.squeeze(metric))
    
        @tullio b_mnl[p,q,Q] := Ppq[p,q,P]*metric[P,Q]
    end

    if(p.gpu)
        if (!p.RI)
            I = cu(I)
        else
            b_mnl = cu(b_mnl)
        end
    end

    return S,T,V,H,I,b_mnl

end


######################################### J_mn^(j) K_mn^(j) #########################################

function computeJKj(C,I,b_mnl,p)

    if(p.gpu)
        if(p.RI)
            J,K = JKj_RI(cu(C),b_mnl,p.nbf,p.nbf5,p.nbfaux)
        else
            J,K = JKj_Full(cu(C),I,p.nbf,p.nbf5)
        end
        return Array(J),Array(K)
    else
        if(p.RI)
            J,K = JKj_RI(C,b_mnl,p.nbf,p.nbf5,p.nbfaux)
        else
            J,K = JKj_Full(C,I,p.nbf,p.nbf5)
        end
        return J,K
    end

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

######################################### J_pq K_pq #########################################

function computeJKH_MO(C,H,I,b_mnl,p)

    if(p.gpu)
        if(p.RI)
            J_MO,K_MO,H_core = JKH_MO_RI(cu(C),cu(H),b_mnl,p.nbf,p.nbf5,p.nbfaux)
        else
            J_MO,K_MO,H_core = JKH_MO_Full(cu(C),cu(H),I,p.nbf,p.nbf5)
        end
        return Array(J_MO),Array(K_MO),Array(H_core)
     else
        if(p.RI)
             J_MO,K_MO,H_core = JKH_MO_RI(C,H,b_mnl,p.nbf,p.nbf5,p.nbfaux)
         else
             J_MO,K_MO,H_core = JKH_MO_Full(C,H,I,p.nbf,p.nbf5)
         end
         return J_MO,K_MO,H_core
      end

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

function JK_HF_Full(D,I,p)

    #denmatj
    @tullio J[m,n] := D[l,s]*I[m,n,s,l]
    @tullio K[m,s] := D[n,l]*I[m,n,s,l]

    return J,K

    end

function iajb_Full_jit(C,I,no1,nalpha,nbf,nbf5)

    Cocc = view(C,:,no1+1:nalpha)
    Ccwo = view(C,:,nalpha+1:nbf)

    @tullio iajb[i,a,j,b] := Cocc[m,i]*Ccwo[n,a]*I[m,n,s,l]*Cocc[s,j]*Ccwo[l,b]

    return iajb

end

end
