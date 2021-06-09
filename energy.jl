module energy
using LinearAlgebra
using Printf

include("minimization.jl")
include("guess.jl")

function compute_energy(S,T,V,H,I,b_mnl,E_nuc,p;fmiug0=nothing,gamma=nothing,hfidr=true,printmode=true)

    Ei,C = eigen(H, S)

    if hfidr
        EHF,Cguess,fmiug0guess = minimization.hfidr(C,H,I,b_mnl,E_nuc,p)
    end

    if gamma==nothing
        gamma = guess.compute_gamma(p)
    end

    C = Cguess
    elag = zeros(p.nbf,p.nbf)
    gamma,n,cj12,ck12 = minimization.occoptr(gamma,true,false,C,H,I,b_mnl,p)

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
        convgdelag,E_old,E_diff,sumdiff_old,itlim,fmiug0,C,elag = minimization.orboptr(C,n,H,I,b_mnl,cj12,ck12,E_old,E_diff,sumdiff_old,i_ext,itlim,fmiug0,E_nuc,p,printmode)
        #t2 = time()

        #occopt
        gamma,n,cj12,ck12 = minimization.occoptr(gamma,false,convgdelag,C,H,I,b_mnl,p)
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
        for i in 1:p.nbf5
            @printf(" %3i    %9.7f  %10.8f",i,2*n[i],elag[i,i])
            println(" ")
        end

        println("----------------")
        println(" Final Energies ")
        println("----------------")

        if hfidr
            @printf("       HF Total Energy = %15.7f\n",E_nuc + EHF)
	end
        @printf("Final NOF Total Energy = %15.7f\n",E_nuc + E_old)
        if hfidr
            @printf("    Correlation Energy = %15.7f\n",E_old-EHF)
	end
        println(" ")
        println(" ")
    end

end

end
