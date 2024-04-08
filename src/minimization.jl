function hfidr(C,H,I,b_mnl,E_nuc,p;printmode=true)

    no1_ori = p.no1
    p.no1 = p.nbeta

    n = zeros((p.nbf5))
    n[1:p.nbeta] .= 1.0
    n[p.nbeta+1:p.nalpha] .= 0.5

    @tullio cj12[i,j] := 2*n[i]*n[j]
    @tullio ck12[i,j] := n[i]*n[j]
    if p.MSpin==0 && p.nsoc>1
        ck12[p.nbeta+1:p.nalpha,p.nbeta+1:p.nalpha] .*= 2
    end

    if(printmode)
        @printf("Hartree-Fock\n")
        @printf("============\n")
        @printf("\n")
        @printf("  %6s  %6s %8s %13s %15s %16s\n","Nitext","Nitint","Eelec","Etot","Ediff","maxdiff")
    end

    E,elag,sumdiff,maxdiff = ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)
    fmiug0 = nothing

    ext = true
    # iteraciones externas
    for i_ext in 1:p.maxitid
        if i_ext==1
            maxlp = 1
        else
            maxlp = p.maxloop
	end

        # iteraciones internas
	E_diff = 0
        for i_int in 1:maxlp
            E_old = E

            if p.scaling
                fmiug = fmiug_scaling(fmiug0,elag,i_ext,p.nzeros,p.nbf,p.noptorb)
            end
            fmiug0, W = eigen(fmiug)
	    @tullio Cnew[i,j] := C[i,k]*W[k,j]
            C = Cnew
            E,elag,sumdiff,maxdiff = ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)

            E_diff = E-E_old
            if(abs(E_diff)<p.thresheid)
                if(printmode)
                    @printf("%6i %6i %14.8f %14.8f %14.8f %14.8f \n",i_ext,maxlp,E,E+E_nuc,E_diff,maxdiff)
                    flush(stdout)
		end
		fmiug0 = diag(elag)
                ext = false
                break
	    end
	end

        if !ext
            break
	end
        if printmode
            @printf("%6i %6i %14.8f %14.8f %14.8f %14.8f \n",i_ext,maxlp,E,E+E_nuc,E_diff,maxdiff)
            flush(stdout)
        end

    end
    # Regresamos no1 a su estado original
    p.no1 = no1_ori

    return E,C,fmiug0

end

function occoptr(gamma,C,H,I,b_mnl,freeze_occ,p)

    if p.ndoc>0 && !freeze_occ
        J_MO,K_MO,H_core = computeJKH_MO(C,H,I,b_mnl,p)
	res = optimize(gamma->calcocce(gamma,J_MO,K_MO,H_core,p),gamma->calcoccg(gamma,J_MO,K_MO,H_core,p),gamma, ConjugateGradient(), Optim.Options(g_abstol = p.threshen), inplace=false)
	gamma = res.minimizer
    end

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
    cj12,ck12 = PNOFi_selector(n,p)

    if p.ndoc>0 && !freeze_occ
        return res.minimum, res.iterations, res.ls_success, gamma, n, cj12, ck12
    else
        return -1, 0, true, gamma, n, cj12, ck12
    end
end

function orboptr(C,n,H,I,b_mnl,cj12,ck12,i_ext,itlim,fmiug0,p,printmode=true)

    i_int = 0
    success_orb = false

    E,elag,sumdiff,maxdiff = ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)

    if maxdiff<p.threshl
        success_orb = true
        return E,C,i_int,success_orb,itlim,fmiug0
    end

    if p.scaling && i_ext>2 && i_ext >= itlim #&& sumdiff > sumdiff_old
        p.nzeros = p.nzeros + 1
        itlim = i_ext + p.itziter
        if p.nzeros>p.nzerosm
            p.nzeros = p.nzerosr
        #if p.nzeros>abs(trunc(Int,log10(maxdiff)))+1
        #    p.nzeros = p.nzerosr#abs(trunc(Int,log10(maxdiff)))
	end
    end
    sumdiff_old = sumdiff

    if fmiug0===nothing && i_ext==1
        maxlp = 1
    else
        maxlp = p.maxloop
    end

    fmiug = zeros(p.noptorb,p.noptorb)
    fk = zeros(30,p.noptorb,p.noptorb)
    bdiis = zeros(31,31)
    cdiis = zeros(31)
    iloop = 0
    idiis = 0

    for i_int in 1:maxlp
        iloop = iloop + 1
        E_old2 = E

        #scaling
        if p.scaling
            fmiug = fmiug_scaling(fmiug0,elag,i_ext,p.nzeros,p.nbf,p.noptorb)
	end
        if p.diis && maxdiff < p.thdiis
            fk,fmiug,idiis,bdiis = fmiug_diis(fk,fmiug,idiis,bdiis,cdiis,maxdiff,p)
        end
        eigval, eigvec = eigen(fmiug)
        fmiug0 = eigval

	@tullio Cnew[i,j] := C[i,k]*eigvec[k,j]
	C = Cnew

        E,elag,sumdiff,maxdiff = ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)

        E_diff2 = E-E_old2

        if abs(E_diff2)<p.threshec || i_int==maxlp
	    break
        end
    end
    return E,C,iloop,success_orb,itlim,fmiug0
end

function orbopt_rotations(gamma,C,H,I,b_mnl,p)

    y = zeros(p.nvar)

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
    cj12,ck12 = PNOFi_selector(n,p)

    res = optimize(Optim.only_fg!((F,G,y)->calcorbeg(F,G,y,n,cj12,ck12,C,H,I,b_mnl,p)), y, ConjugateGradient(), Optim.Options(iterations = p.maxloop), inplace=false)

    E = res.minimum
    y = res.minimizer

    C = rotate_orbital(y,C,p)
    nit = res.iterations
    success = res.g_converged

    return E,C,nit,success
end

function comb(gamma,C,H,I,b_mnl,p)

    x = zeros(p.nvar+p.nv)
    x[p.nvar+1:end] = gamma

    res = optimize(Optim.only_fg!((F,G,x)->calccombeg(F,G,x,C,H,I,b_mnl,p)), x, ConjugateGradient(), Optim.Options(iterations = p.maxloop), inplace=false)

    E = res.minimum
    x = res.minimizer
    y = x[1:p.nvar]
    gamma = x[p.nvar+1:end]
    C = rotate_orbital(y,C,p)

    n,dR = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)

    return E,C,gamma,n,res.iterations,res.g_converged
end
