function compute_geom_gradients(bset,n,C,cj12,ck12,elag,p)

    C_nbf5 = @view C[1:end,1:p.nbf5]

    @tullio RDM1[m,l] := 2 * n[p] * C_nbf5[m,p] * C_nbf5[l,p]
    @tullio lag[m,l] := 2 * C[m,q] * elag[q,p] *C[l,p]

    grad = compute_grad_nuc(bset, p)

    for i in 1:p.natoms
	grad_i = @view grad[i, 1:end]

	dS = ∇overlap(bset, i)
	@tullio grad_i[k] += -lag[m,n] * dS[m,n,k]

	dT = ∇kinetic(bset, i)
	@tullio grad_i[k] += RDM1[m,n] * dT[m,n,k]

	dV = ∇nuclear(bset, i)
	@tullio grad_i[k] += RDM1[m,n] * dV[m,n,k]
    end

    cj12[diagind(cj12)] .= 0
    ck12[diagind(ck12)] .= 0

    n_beta = @view n[1:p.nbeta]
    n_bf5 = @view n[p.nalpha+1:p.nbf5]

    C_beta = @view C[1:end, 1:p.nbeta]
    C_bf5 = @view C[1:end, p.nalpha+1:p.nbf5]

    if(p.RI == false)
        @tullio RDM2[m,n,s,l] :=  cj12[p,q] * C_nbf5[m,p] * C_nbf5[n,p] * C_nbf5[s,q] * C_nbf5[l,q]
        @tullio RDM2[m,n,s,l] +=  n_beta[p] * C_beta[m,p] * C_beta[n,p] * C_beta[s,p] * C_beta[l,p]
        @tullio RDM2[m,n,s,l] +=  n_bf5[p]  * C_bf5[m,p] * C_bf5[n,p] * C_bf5[s,p] * C_bf5[l,p]
        @tullio RDM2[m,n,s,l] += -ck12[p,q] * C_nbf5[m,p] * C_nbf5[l,p] * C_nbf5[s,q] * C_nbf5[n,q]

        for i in 1:p.natoms
            grad_i = @view grad[i, 1:end]

            dERI = ∇ERI_2e4c(bset, i)
            @tullio grad_i[k] += RDM2[m,n,s,l] * dERI[m,n,s,l,k]
        end
    else

        aux = nothing
        try
            aux = BasisSet(bset.name * "-jkfit", bset.atoms)
        catch
            aux = BasisSet("def2-universal-jkfit", bset.atoms)
        end

        G = ERI_2e2c(aux)
        G = (G .+ G') ./ 2
        I = ERI_2e3c(bset, aux)
        G = inv(G)

        nbf = size(I)[1]
        naux = size(G)[1]

        @tullio IG[m,n,Q] := I[m,n,P] * G[P,Q]
        I = nothing
        G = nothing

        @tullio I2[p,q,Q] := C_nbf5[m,p] * C_nbf5[n,q] * IG[m,n,Q]

        @tullio bk[p,q,Q] := -ck12[p,q] * I2[p,q,Q]
        @tullio d[l,s,Q] := C_nbf5[l,p] * C_nbf5[s,q] * bk[p,q,Q]
	bk = nothing

        @tullio bj[q,Q] := cj12[p,q] * I2[p,p,Q]
        @tullio d[l,s,Q] += C_nbf5[l,q] * C_nbf5[s,q] * bj[q,Q]
	bj = nothing

        I_beta = @view I2[1:p.nbeta,1:p.nbeta,1:end]
        @tullio bn_beta[p,Q] := n_beta[p] * I_beta[p,p,Q]
        @tullio d[l,s,Q] += C_beta[l,q] * C_beta[s,q] * bn_beta[q,Q]
	bn_beta = nothing

        I_bf5 = @view I2[p.nalpha+1:end,p.nalpha+1:end,1:end]
        @tullio bn_bf5[p,Q] := n_bf5[p] * I_bf5[p,p,Q]
        @tullio d[l,s,Q] += C_bf5[l,q] * C_bf5[s,q] * bn_bf5[q,Q]
	bn_bf5 = nothing

        I2 = nothing

	@tullio m[Q,R] := d[l,s,Q] * IG[s,l,R]


        for i in 1:p.natoms
            geri3c = ∇ERI_2e3c(bset, aux, i)
            geri2c = ∇ERI_2e2c(aux, i)

            grad_i = @view grad[i, 1:end]

            @tullio grad_i[k] += 2 * d[s,l,Q] * geri3c[s,l,Q,k]
            @tullio grad_i[k] += - m[Q,R] * geri2c[Q,R,k] 

        end
    end

    return grad
end
