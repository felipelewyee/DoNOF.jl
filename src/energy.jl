export energy

function energy(bset,p;C=nothing,fmiug0=nothing,n=nothing,do_hfidr=true,do_nofmp2=false,printmode=true,nofmp2strategy="numerical",tolnofmp2=1e-10,do_ekt=false,do_mulliken_pop=false,do_lowdin_pop=true,do_m_diagnostic=true,do_mbpt=false,freeze_occ=false,do_translate_to_donofsw=false,do_erpa=false)

    t0 = time()
    println("Computing Integrals")
    flush(stdout)

    if p.gpu
        S,T,V,H,I,b_mnl = compute_integrals(bset,p,true)
    else
        S,T,V,H,I,b_mnl = compute_integrals(bset,p)
    end

    t1 = time()
    @printf("Elapsed time: %7.2f Seconds\n\n", t1-t0)
    flush(stdout)

    if(printmode)
        println("Number of basis functions                   (NBF)    = ",p.nbf)
        if(p.RI)
            println("Number of auxiliary basis functions         (NBFAUX) = ",p.nbfaux)
	end
        println("Inactive Doubly occupied orbitals up to     (NO1)    = ",p.no1)
        println("No. considered Strongly Doubly occupied MOs (NDOC)   = ",p.ndoc)
        println("No. considered Strongly Singly occupied MOs (NSOC)   = ",p.nsoc)
        println("NO. of Weakly occ. per St. Doubly occ.  MOs (NCWO)   = ",p.ncwo)
        println("Dimension of the Nat. Orb. subspace         (NBF5)   = ",p.nbf5)
        println("No. of electrons                                     = ",p.ne)
        println("    No. of alpha electrons                           = ",p.nalpha)
        println("    No. of beta electrons                            = ",p.nbeta)
        println("Multiplicity                                         = ",p.mul)
        println("")
    end
    flush(stdout)

    println("Geometry")
    println("========")
    for i in 1:bset.natoms
	    @printf("%2s %15.7f %15.7f %15.7f\n",Z_to_symbol(bset.atoms[i].Z),bset.atoms[i].xyz[1],bset.atoms[i].xyz[2],bset.atoms[i].xyz[3])
    end
    println("")

    E_nuc = compute_E_nuc(bset,p)

    if isnothing(C)
        Ei,Cguess = eigen(H, S)
    else
        Cguess = C
    end
    Cguess = check_ortho(Cguess,S,p)

    if do_hfidr
        EHF,Cguess,fmiug0guess = hfidr(Cguess,H,I,b_mnl,E_nuc,p)
    end

    if isnothing(C)
        C = Cguess
    end
    C = check_ortho(C,S,p)

    if isnothing(n)
        gamma = guess_gamma_trigonometric(p)
        n,dn_dgamma = ocupacion_trigonometric(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin)
        if p.occ_method == "Softmax"
	    p.nv = p.nbf5 - p.no1 - p.nsoc
	    gamma = n_to_gammas_softmax(n,p)
        end
    else
        if p.occ_method == "Trigonometric"
            gamma = n_to_gammas_trigonometric(n,p)
        elseif p.occ_method == "Softmax"
	    p.nv = p.nbf5 - p.no1 - p.nsoc
	    gamma = n_to_gammas_softmax(n,p)
	end
    end

    elag = zeros(p.nbf,p.nbf)
    E_occ,nit_occ,success_occ,gamma,n,cj12,ck12 = occoptr(gamma,C,H,I,b_mnl,freeze_occ,p)

    itlim = 1
    E = 9999#EHF
    E_old = 9999#EHF
    E_diff = 9999
    sumdiff_old = 0

    if printmode 
        println(" ")
        @printf("PNOF%i Calculation\n",p.ipnof)
        println("==================")
        println(" ")
    end

    @printf("  %6s  %7s %5s   %7s %5s      %8s %13s %15s %12s    %8s  %8s\n","Nitext","Nit_orb","Time","Nit_occ","Time","Eelec","Etot","Ediff","maxdiff","Grad_orb","Grad_occ")

    for i_ext in 1:p.maxit
	ta1 = time()
        if p.orb_method == "ID"
            E_orb,C,nit_orb,success_orb,itlim,fmiug0 = orboptr(C,n,H,I,b_mnl,cj12,ck12,i_ext,itlim,fmiug0,p,printmode)
        elseif p.orb_method == "Rotations"
            E_orb,C,nit_orb,success_orb = orbopt_rotations(gamma,C,H,I,b_mnl,p)
        end
	ta2 = time()
    
        E_occ,nit_occ,success_occ,gamma,n,cj12,ck12 = occoptr(gamma,C,H,I,b_mnl,freeze_occ,p)
	if(p.occ_method=="Softmax")
            C,gamma = order_occupations_softmax(C,gamma,H,I,b_mnl,p)
	end
	ta3 = time()

        E = E_orb
        E_diff = E-E_old
        E_old = E

        Etmp,elag,sumdiff,maxdiff = ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)

        y = zeros(p.nvar)
        n,dn_dgamma = ocupacion(gamma,p.no1,p.ndoc,p.nalpha,p.nv,p.nbf5,p.ndns,p.ncwo,p.HighSpin,p.occ_method)
        cj12,ck12 = PNOFi_selector(n,p)
        grad_orb = calcorbg(y,n,cj12,ck12,C,H,I,b_mnl,p)
        J_MO,K_MO,H_core = computeJKH_MO(C,H,I,b_mnl,p)
        grad_occ = calcoccg(gamma,J_MO,K_MO,H_core,p)

        M = M_diagnostic(p,n,get_value=true)
        @printf("%6i %7i %10.1e %4i %10.1e %14.8f %14.8f %14.8f %10.6f   %4.1e   %4.1e %4.2f\n",i_ext,nit_orb,ta2-ta1,nit_occ,ta3-ta2,E,E+E_nuc,E_diff,maxdiff,norm(grad_orb),norm(grad_occ),M)

        if isnothing(fmiug0)
            save(p.title*".jld", "E", Etmp, "C", C,"n",n)
        else
            save(p.title*".jld", "E", Etmp, "C", C,"n",n,"fmiug0",fmiug0)
        end
        flush(stdout)

        if(norm(grad_orb) < p.threshgorb && norm(grad_occ) < p.threshgocc)
            break
        end

    end

    C,n,elag = order_subspaces(C,n,elag,H,I,b_mnl,p)
    if isnothing(fmiug0)
        save(p.title*".jld", "C", C,"n",n)
    else
        save(p.title*".jld", "C", C,"n",n,"fmiug0",fmiug0)
    end

    if printmode
        println(" ")
        println("RESULTS OF THE OCCUPATION OPTIMIZATION")
        println("========================================")
        for i in 1:p.nbeta
            @printf(" %3i    %9.7f  %10.8f\n",i,2*n[i],elag[i,i])
        end
        for i in p.nbeta+1:p.nalpha
	    if !p.HighSpin
                @printf(" %3i    %9.7f  %10.8f\n",i,2*n[i],elag[i,i])
	    else
                @printf(" %3i    %9.7f  %10.8f\n",i,n[i],elag[i,i])
	    end
        end
        for i in p.nalpha+1:p.nbf5
            @printf(" %3i    %9.7f  %10.8f\n",i,2*n[i],elag[i,i])
        end

        println(" ")
        println("----------------")
        println(" Final Energies ")
        println("----------------")

        if do_hfidr
            @printf("       HF Total Energy = %15.7f\n",E_nuc + EHF)
	end
        @printf("Final NOF Total Energy = %15.7f\n",E_nuc + E)
        if do_hfidr
            @printf("    Correlation Energy = %15.7f\n",E-EHF)
	end
        println(" ")
        println(" ")
    end

    E_t = E_nuc + E

    fchk(p.title,p,bset,"Energy",E_t,elag,n,C)

    t2 = time()
    @printf("Elapsed time: %7.2f Seconds\n", t2-t1)
    flush(stdout)

    if do_nofmp2
        nofmp2(n,C,H,I,b_mnl,E_nuc,p,nofmp2strategy,tolnofmp2)
    end

    if(do_ekt)
        ext_koopmans(p,elag,n)
    end

    if(do_mulliken_pop)
        mulliken_pop(bset,p,n,C,S)
    end

    if(do_lowdin_pop)
        lowdin_pop(bset,p,n,C,S)
    end

    if(do_m_diagnostic)
        M_diagnostic(p,n)
    end

    if(do_mbpt)
        mbpt(n,C,H,I,b_mnl,E_nuc,E,p)
    end

    if(do_erpa)
        erpa(n,C,H,I,b_mnl,E_nuc,E,p)
    end

    if(do_translate_to_donofsw)
        write_to_DoNOFsw(p,bset,n,C,diag(elag),fmiug0,10,E_nuc + E)
    end

    return E_nuc + E,C,n,fmiug0

end
