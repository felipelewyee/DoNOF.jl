module integrals

using Tullio
using LoopVectorization

function JKj_Full(C,I,b_mnl,p)

    Cnbf5 = view(C,:,1:p.nbf5)
    
    @tullio D[i,m,n] := Cnbf5[m,i]*Cnbf5[n,i]
    @tullio J[i,m,n] := D[i,s,l]*I[m,n,s,l]
    @tullio K[i,m,s] := D[i,n,l]*I[m,n,s,l]

    return J,K

end

function JKH_MO_Full(C,H,I,p)

    #denmatj

    Cnbf5 = view(C,:,1:p.nbf5)

    @tullio D[i,m,n] := Cnbf5[m,i]*Cnbf5[n,i]
    #D = np.einsum('mi,ni->imn', C[:,0:p.nbf5], C[:,0:p.nbf5],optimize=True)
    #QJMATm
    @tullio J[j,m,n] := D[j,s,l]*I[m,n,s,l]
    #J = np.tensordot(D, I, axes=([1,2],[2,3]))
    @tullio J_MO[i,j] := D[i,m,n]*J[j,m,n]
    #J_MO = np.tensordot(J, D,axes=((1,2),(1,2)))
    #QKMATm
    @tullio K[j,m,s] := D[j,n,l]*I[m,n,s,l]
    #K = np.tensordot(D, I, axes=([1,2],[1,3]))
    @tullio K_MO[i,j] := D[i,m,s]*K[j,m,s]
    #K_MO = np.tensordot(K, D, axes=([1,2],[1,2]))
    #QHMATm
    @tullio H_core[i] := D[i,m,n]*H[m,n]
    #H_core = np.tensordot(D, H, axes=([1,2],[0,1]))

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

end
