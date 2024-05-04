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

    n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
    cj12,ck12 = PNOFi_selector(n,p)

    if p.orb_method=="Rotations"

        y = zeros(p.nvar)
        res = optimize(Optim.only_fg!((F,G,y)->calcorbeg(F,G,y,n,cj12,ck12,C,H,I,b_mnl,p)), y, ConjugateGradient(), Optim.Options(iterations = p.maxloop), inplace=false)
    
        E = res.minimum
        y = res.minimizer
    
        C = rotate_orbital(y,C,p)
        nit = res.iterations
        success = res.g_converged

    elseif p.orb_method=="ADAM"  
        E,C,nit,success = experimental_minimize_rotations(n,cj12,ck12,C,H,I,b_mnl,p)
    end

    return E,C,nit,success
end

#### Experimental ####
# ADAM
function experimental_minimize_rotations(n,cj12,ck12,C,H,I_AO,b_mnl,p)


    elag,Hmat = compute_Lagrange2(C,n,H,I_AO,b_mnl,cj12,ck12,p)
    E = computeE_elec(Hmat,n,elag,p)

    alpha = p.alpha
    beta1 = 0.7
    beta2 = 0.999

    y = zeros(p.nvar)
    m = 0.0 .* y
    v = 0.0 .* y
    vhat_max = 0.0 .* y

    improved = false
    success = false
    best_E, best_C = E, C
    nit = 0
    for i in 1:p.maxloop
        nit = nit + 1

        grad = 4*elag - 4*elag'
        grads = zeros(p.nvar)
        nn = 1
        for i in 1:p.nbf5
            for j in i+1:p.nbf
                grads[nn] = grad[i,j]
                nn += 1
            end
        end

	#if(i>1)
	#  println(maximum(abs.(grads)), " ", norm(grads))
        #end

        if norm(grads) < p.threshgorb && improved
            success = true
            break
        end

        m = beta1 .* m + (1 - beta1) .* grads
	v = beta2 .* v + (1 - beta2) .* (grads .^ 2)
	mhat = m ./ (1.0 - beta1^i)
	vhat = v ./ (1.0 - beta2^i)
	vhat_max = max.(vhat_max, vhat)
	y = - alpha * mhat ./ (sqrt.(vhat_max .+ 10^-8)) #AMSgrad
        C = rotate_orbital(y,C,p)

        elag,Hmat = compute_Lagrange2(C,n,H,I_AO,b_mnl,cj12,ck12,p)
        E = computeE_elec(Hmat,n,elag,p)
	#println(i," ",E," ", E < best_E)
        if E < best_E
            best_C = C
            best_E = E
	    improved = true
        end

    end

    if !improved
        p.alpha = p.alpha/10
        p.maxloop = p.maxloop + 30
	#println("      alpha ",p.alpha)
    end

    return best_E,best_C,nit,success
end

# RMSprop
#function experimental_minimize_rotations(n,cj12,ck12,C,H,I_AO,b_mnl,p)
#
#    best_C = C
#
#    y = zeros(p.nvar)
#    nit = 0
#    elag,Hmat = compute_Lagrange2(C,n,H,I_AO,b_mnl,cj12,ck12,p)
#    E = computeE_elec(Hmat,n,elag,p)
#    E_old = E
#    best_E = E
#    success = false
#    sg = 0.0 .* y
#    alpha = 0.0 .* y
#    eta = 0.0004
#    rho = 0.999
#    for i in 1:p.maxloop
#        nit = nit + 1
#        C_new = rotate_orbital(y,C,p)
#
#        elag,Hmat = compute_Lagrange2(C_new,n,H,I_AO,b_mnl,cj12,ck12,p)
#
#        E = computeE_elec(Hmat,n,elag,p)
#        if E <= best_E
#            best_C = C_new
#            best_E = E
#        end
#
#        println(E," ", E <= best_E)
#        E_old = E
#        C = C_new
#
#        grad = 4*elag - 4*elag'
#        grads = zeros(p.nvar)
#        nn = 1
#        for i in 1:p.nbf5
#            for j in i+1:p.nbf
#                grads[nn] = grad[i,j]
#                nn += 1
#            end
#        end
#        if norm(grads) < p.threshgorb
#            success = true
#            break
#        end
#
#        sg = rho .* sg + (1-rho) .* (grads .^ 2)
#	alpha = eta ./ sqrt.(sg .+ 10^-8)
#	y = -alpha .* grads
#    end
#
#    elag,Hmat = compute_Lagrange2(C,n,H,I_AO,b_mnl,cj12,ck12,p)
#    E = computeE_elec(Hmat,n,elag,p)
#    return best_E,best_C,nit,success
#end

# Momentum Gradient Descent
#function experimental_minimize_rotations(n,cj12,ck12,C,H,I_AO,b_mnl,p)
#
#    best_C = C
#
#    y = zeros(p.nvar)
#    nit = 0
#    elag,Hmat = compute_Lagrange2(C,n,H,I_AO,b_mnl,cj12,ck12,p)
#    E = computeE_elec(Hmat,n,elag,p)
#    E_old = E
#    best_E = E
#    success = false
#    step = 0 .* y
#    momentum = 0.97
#    alpha = 1.0-momentum
#    for i in 1:p.maxloop
#	nit = nit + 1
#        C_new = rotate_orbital(y,C,p)
#
#        elag,Hmat = compute_Lagrange2(C_new,n,H,I_AO,b_mnl,cj12,ck12,p)
#        
#        E = computeE_elec(Hmat,n,elag,p)
#	if E <= best_E
#	    best_C = C_new
#	    best_E = E
#	end
#
#        println(E," ",alpha," ", E <= best_E)
#        E_old = E
#  	C = C_new
#
#        grad = 4*elag - 4*elag'
#        grads = zeros(p.nvar)
#        nn = 1
#        for i in 1:p.nbf5
#            for j in i+1:p.nbf
#                grads[nn] = grad[i,j]
#                nn += 1
#            end
#        end
#	if norm(grads) < p.threshgorb
#	    success = true
#	    break
#	end
#	step = -alpha*grads + momentum*step
#        y = step
#    end
#
#    elag,Hmat = compute_Lagrange2(C,n,H,I_AO,b_mnl,cj12,ck12,p)
#    E = computeE_elec(Hmat,n,elag,p)
#    return best_E,best_C,nit,success
#end

######################

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
