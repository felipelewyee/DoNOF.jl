function compute_energy(wfn,mol,p;C=nothing,fmiug0=nothing,gamma=nothing,do_hfidr=true,do_nofmp2=false,printmode=true)

    S,T,V,H,I,b_mnl = compute_integrals(wfn,mol,p)

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
        println("Multiplicity                                         = ",p.mul)
        println("")
    end

    E_nuc = mol.nuclear_repulsion_energy()

    Cguess = C
    if isnothing(C)
        Ei,Cguess = eigen(H, S)
    end

    if do_hfidr
        EHF,Cguess,fmiug0guess = hfidr(Cguess,H,I,b_mnl,E_nuc,p)
    end

    if isnothing(C)
        C = Cguess
    end

    if isnothing(gamma)
        gamma = compute_gamma(p)
    end

    C = Cguess
    elag = zeros(p.nbf,p.nbf)
    gamma,n,cj12,ck12 = occoptr(gamma,true,false,C,H,I,b_mnl,p)

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
    for i_ext in 1:p.maxit
        #t1 = time()
        #orboptr
        convgdelag,E_old,E_diff,sumdiff_old,itlim,fmiug0,C,elag = orboptr(C,n,H,I,b_mnl,cj12,ck12,E_old,E_diff,sumdiff_old,i_ext,itlim,fmiug0,E_nuc,p,printmode)
        #t2 = time()

        #occopt
        gamma,n,cj12,ck12 = occoptr(gamma,false,convgdelag,C,H,I,b_mnl,p)
        #t3 = time()

        if convgdelag
            break
	end
        #print(t2-t1,t3-t2)
    end

    if printmode
        println(" ")
        println("RESULTS OF THE OCCUPATION OPTIMIZATION")
        println("========================================")
        for i in 1:p.nbeta
            @printf(" %3i    %9.7f  %10.8f\n",i,2*n[i],elag[i,i])
        end
        for i in p.nbeta+1:p.nalpha
	    if p.MSpin==0
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
        nofmp2(n,C,H,I,b_mnl,E_nuc,p)
    end

    return E_nuc + E_old,C,gamma,fmiug0

end
