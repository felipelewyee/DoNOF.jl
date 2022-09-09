function energy(bset,p;C=nothing,fmiug0=nothing,gamma=nothing,do_hfidr=true,do_nofmp2=false,printmode=true,nofmp2strategy="numerical",tolnofmp2=1e-10,do_ekt=false,do_mulliken_pop=false,do_lowdin_pop=false,do_m_diagnostic=false,freeze_occ=false)

    t1 = time()

    S,T,V,H,I,b_mnl = compute_integrals(bset,p)

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

    println("Geometry")
    println("========")
    for i in 1:bset.natoms
	    @printf("%2s %15.7f %15.7f %15.7f\n",Z_to_symbol(bset.atoms[i].Z),bset.atoms[i].xyz[1],bset.atoms[i].xyz[2],bset.atoms[i].xyz[3])
    end
    println("")

    E_nuc = compute_E_nuc(bset,p)

    Cguess = C
    if isnothing(C)
        Ei,Cguess = eigen(H, S)
    end
    Cguess = check_ortho(Cguess,S,p)

    if do_hfidr
        EHF,Cguess,fmiug0guess = hfidr(Cguess,H,I,b_mnl,E_nuc,p)
    end

    if isnothing(C)
        C = Cguess
    end
    Cguess = check_ortho(C,S,p)

    if isnothing(gamma)
        gamma = compute_gamma(p)
    end

    C = Cguess
    elag = zeros(p.nbf,p.nbf)
    gamma,n,cj12,ck12 = occoptr(gamma,C,H,I,b_mnl,freeze_occ,p)

    iloop = 0
    itlim = 1
    E_old = 9999#EHF
    E_diff = 9999
    sumdiff_old = 0

    if printmode 
        println(" ")
        @printf("PNOF%i Calculation\n",p.ipnof)
        println("==================")
        println(" ")
        @printf("  %6s  %6s %8s %13s %15s %16s\n","Nitext","Nitint","Eelec","Etot","Ediff","maxdiff")
    end

    if p.method == "ID"

        for i_ext in 1:p.maxit
	    ta1 = time()
            convgdelag,E_old,E_diff,sumdiff_old,itlim,fmiug0,C,elag = orboptr(C,n,H,I,b_mnl,cj12,ck12,E_old,E_diff,sumdiff_old,i_ext,itlim,fmiug0,E_nuc,p,printmode)
	    ta2 = time()
    
            gamma,n,cj12,ck12 = occoptr(gamma,C,H,I,b_mnl,freeze_occ,p)
	    ta3 = time()
	    @printf("Orb: %6.2e Occ: %6.2e\n", ta2-ta1, ta3-ta2)

            if convgdelag
                break
    	    end
            save(p.title*".jld","C", C,"gamma",gamma,"fmiug0",fmiug0)
        end
    end

    if p.method == "Rotations"

	E_old = 0
        for i_ext in 1:p.maxit
	    ta1 = time()
            E,C,nit_orb,success_orb = orbopt_rotations(gamma,C,H,I,b_mnl,p)
	    ta2 = time()

            gamma,n,cj12,ck12 = occoptr(gamma,C,H,I,b_mnl,freeze_occ,p)
	    ta3 = time()
            Etmp,elag,sumdiff,maxdiff = ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)
	    ta4 = time()

	    @printf("%6i %6i %14.8f %14.8f %14.8f %14.8f\n",i_ext,nit_orb,E,E+E_nuc,E-E_old,maxdiff)
	    @printf("Orb: %6.2e Occ: %6.2e Lag: %6.2e\n", ta2-ta1, ta3-ta2, ta4-ta3)

	    if(abs(E-E_old) < 1e-4)
	        break
	    end
            E_old = E
            save(p.title*".jld","C", C,"gamma",gamma,"fmiug0",fmiug0)
	    flush(stdout)
        end

    end

    save(p.title*".jld","C", C,"gamma",gamma,"fmiug0",fmiug0)

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
        @printf("Final NOF Total Energy = %15.7f\n",E_nuc + E_old)
        if do_hfidr
            @printf("    Correlation Energy = %15.7f\n",E_old-EHF)
	end
        println(" ")
        println(" ")
    end

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

    t2 = time()
    @printf("Elapsed time: %7.2f\n Seconds", t2-t1)

    return E_nuc + E_old,C,gamma,fmiug0

end

function brute_force_energy(bset,p;intents=5,C=nothing,gamma=nothing,fmiug0=nothing,do_hfidr=true,RI_last=false,gpu_last=false,do_ekt=false,do_mulliken_pop=false,do_lowdin_pop=false,do_m_diagnostic=false)

    E,C,gamma,fmiug0 = energy(bset,p,C=C,gamma=gamma,fmiug0=fmiug0,do_hfidr=do_hfidr)
    E_min = E
    C_min = C
    gamma_min = gamma
    fmiug0_min = fmiug0

    for i in 1:intents
        autozeros(p)

        E,C,gamma,fmiug0 = energy(bset,p,C=C,gamma=gamma,fmiug0=nothing,do_hfidr=false)
        if(E<E_min)
            E_min = E
            C_min = C
            gamma_min = gamma
            fmiug0_min = fmiug0
        end
    end

    p.RI = RI_last
    p.gpu = gpu_last
    E,C,gamma,fmiug0 = energy(bset,p,C=C_min,gamma=gamma_min,fmiug0=fmiug0_min,do_hfidr=false,do_ekt=do_ekt,do_mulliken_pop=do_mulliken_pop,do_lowdin_pop=do_lowdin_pop,do_m_diagnostic=do_m_diagnostic)

    @printf("Best Total NOF Energy %15.7f\n",E)

    return E,C,gamma,fmiug0

end
