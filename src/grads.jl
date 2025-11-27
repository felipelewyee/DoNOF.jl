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

    @tullio RDM2[m,n,s,l] :=  cj12[p,q] * C_nbf5[m,p] * C_nbf5[n,p] * C_nbf5[s,q] * C_nbf5[l,q]
    @tullio RDM2[m,n,s,l] +=  n_beta[p] * C_beta[m,p] * C_beta[n,p] * C_beta[s,p] * C_beta[l,p]
    @tullio RDM2[m,n,s,l] +=  n_bf5[p]  * C_bf5[m,p] * C_bf5[n,p] * C_bf5[s,p] * C_bf5[l,p]
    @tullio RDM2[m,n,s,l] += -ck12[p,q] * C_nbf5[m,p] * C_nbf5[l,p] * C_nbf5[s,q] * C_nbf5[n,q]

    for i in 1:p.natoms
	grad_i = @view grad[i, 1:end]

        dERI = ∇ERI_2e4c(bset, i)
	@tullio grad_i[k] += RDM2[m,n,s,l] * dERI[m,n,s,l,k]

    end

    return grad
end
