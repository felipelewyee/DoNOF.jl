function read_C(;title = "donof")

    C = load(title*".jld")["C"]
    return C

end

function read_n(;title = "donof")

    n = load(title*".jld")["n"]
    return n

end

function read_fmiug0(;title = "donof")
   
    fmiug0 = load(title*".jld")["fmiug0"]
    return fmiug0

end

function read_all(;title = "donof")
    C = read_C(title=title)
    n = read_n(title=title)
    fmiug0 = read_fmiug0(title=title)

    return C,n,fmiug0

end

function increment_C(old_C,old_ncwo,p)

    C = zeros(p.nbf,p.nbf)
    C[1:p.nbf,1:p.no1+p.ndns] = old_C[1:p.nbf,1:p.no1+p.ndns]
    for i in 1:p.ndoc
        old_ll = p.no1 + p.ndns + old_ncwo*(p.ndoc - i) + 1
        old_ul = p.no1 + p.ndns + old_ncwo*(p.ndoc - i + 1)
        ll = p.no1 + p.ndns + p.ncwo*(p.ndoc - i) + 1
        ul = p.no1 + p.ndns + p.ncwo*(p.ndoc - i + 1)
	dif_ncwo = p.ncwo - old_ncwo
        C[1:p.nbf,ll:ll+old_ncwo-1] = old_C[1:p.nbf,old_ll:old_ul]
	C[1:p.nbf,ll+old_ncwo:ul] = old_C[1:p.nbf,p.no1+p.ndns+old_ncwo*p.ndoc+dif_ncwo*(i-1)+1:p.no1+p.ndns+old_ncwo*p.ndoc+dif_ncwo*i]
    end
    C[1:p.nbf,p.no1+p.ndns+p.ndoc*p.ncwo+1:p.nbf] = old_C[1:p.nbf,p.no1+p.ndns+p.ndoc*p.ncwo+1:p.nbf]

    return C

end

function increment_gamma(old_gamma,old_ncwo,p)
    
    gamma = ones(p.nv)*pi/2
    for i in 1:p.ndoc
        gamma[i] = old_gamma[i]
	old_ll = p.ndoc + (old_ncwo-1)*(i-1) + 1
	old_ul = p.ndoc + (old_ncwo-1)*i
	ll = p.ndoc + (p.ncwo-1)*(i - 1) + 1
	ul = p.ndoc + (p.ncwo-1)*i
	gamma[ll:ll+(old_ncwo-1)-1] = old_gamma[old_ll:old_ul]
    end
    
    return gamma
end

function order_subspaces(old_C,old_n,elag,H,I,b_mnl,p)

    C = zeros(p.nbf,p.nbf)
    n = zeros(p.nbf5)
    
    #Sort no1 orbitals
    elag_diag = diag(elag)[1:p.no1]
    sort_idx = sortperm(elag_diag)
    for i in 1:p.no1
	i_old  = sort_idx[i]
        C[1:p.nbf,i] = old_C[1:p.nbf,i_old] 
    end

    #Sort ndoc subspaces
    elag_diag = diag(elag)[p.no1+1:p.ndoc]
    sort_idx = sortperm(elag_diag)
    C[1:p.nbf,1:p.no1] = old_C[1:p.nbf,1:p.no1]
    n[1:p.no1] = old_n[1:p.no1]
    for i in 1:p.ndoc
	i_old  = sort_idx[i]
        ll = p.no1 + p.ndns + p.ncwo*(p.ndoc - i) + 1
        ul = p.no1 + p.ndns + p.ncwo*(p.ndoc - i + 1)
        ll_old = p.no1 + p.ndns + p.ncwo*(p.ndoc - i_old) + 1
        ul_old = p.no1 + p.ndns + p.ncwo*(p.ndoc - i_old + 1)

        C[1:p.nbf,p.no1+i] = old_C[1:p.nbf,p.no1+i_old]
        C[1:p.nbf,ll:ul] = old_C[1:p.nbf,ll_old:ul_old]
        n[p.no1+i] = old_n[p.no1+i_old]
        n[ll:ul] = old_n[ll_old:ul_old]
    end

    #Sort nsoc orbitals
    elag_diag = diag(elag)[p.no1+p.ndoc+1:p.no1+p.ndns]
    sort_idx = sortperm(elag_diag)
    for i in 1:p.nsoc
	i_old  = sort_idx[i]
        C[1:p.nbf,p.no1+p.ndoc+i] = old_C[1:p.nbf,p.no1+p.ndoc+i_old] 
        n[p.no1+p.ndoc+i] = old_n[p.no1+p.ndoc+i_old] 
    end

    C[1:p.nbf,p.nbf5+1:p.nbf] = old_C[1:p.nbf,p.nbf5+1:p.nbf]

    cj12,ck12 = PNOFi_selector(n,p)
    Etmp,elag,sumdiff,maxdiff = ENERGY1r(C,n,H,I,b_mnl,cj12,ck12,p)

    return C,n,elag

end

function write_to_DoNOFsw(p,bset,n,C,elag,fmiug0,it,E)

    Cnew = zeros(p.nbf,p.nbf)

    if(p.spherical)
        println("Warning: Currently DoNOF only support cartesian functions.")
    else
        i = 0
        for basis in bset.basis
            l = basis.l
            ori = trunc(Int,round((l+1)*(l+2)/2))
            if l==2
                Cnew[i+1,1:end] = C[i+1,1:end]
                Cnew[i+4,1:end] = C[i+2,1:end]
                Cnew[i+5,1:end] = C[i+3,1:end]
                Cnew[i+2,1:end] = C[i+4,1:end]
                Cnew[i+6,1:end] = C[i+5,1:end]
                Cnew[i+3,1:end] = C[i+6,1:end]
            elseif l==3
                Cnew[i+1,1:end] = C[i+1,1:end]
                Cnew[i+4,1:end] = C[i+2,1:end]
                Cnew[i+5,1:end] = C[i+3,1:end]
                Cnew[i+6,1:end] = C[i+4,1:end]
                Cnew[i+10,1:end] = C[i+5,1:end]
                Cnew[i+8,1:end] = C[i+6,1:end]
                Cnew[i+2,1:end] = C[i+7,1:end]
                Cnew[i+7,1:end] = C[i+8,1:end]
                Cnew[i+9,1:end] = C[i+9,1:end]
                Cnew[i+2,1:end] = C[i+10,1:end]
            else
                Cnew[i+1:i+ori,1:end] = C[i+1:i+ori,1:end]
            end
	    i += ori
        end
    end

    f = open(p.title*".gcf","w")
    for i in 1:p.nbf5
        @printf(f,"%6d %30.16f\n",i,n[i])
    end
    for i in p.nbf5+1:p.nbf
        @printf(f,"%6d %30.16f\n",i,0.0)
    end
    sumsl = p.nbeta
    for i in 1:p.nbeta
        sumsl -= n[i]
    end
    @printf(f,"%30.16f\n",sumsl)
    for i in 1:p.nbf
        for j in 1:p.nbf
            @printf(f,"%6d %30.16f\n",j,Cnew[j,i])
        end
    end
    for i in 1:p.nbf
        @printf(f,"%6d %30.16f\n",i,elag[i])
    end
    for i in 1:p.nbf
        @printf(f,"%6d %30.16f\n",i,fmiug0[i])
    end
    @printf(f,"%6d %30.16f\n",0,E)
    @printf(f,"%6d %6d %6d %6d %6d %6d\n",p.no1,p.ndoc,p.nsoc,p.ncwo,p.nac,p.no0)
    for i in 1:bset.natoms
	@printf(f,"%6d%30.16f%30.16f%30.16f\n",bset.atoms[i].Z,bset.atoms[i].xyz[1],bset.atoms[i].xyz[2],bset.atoms[i].xyz[3])
    end

    close(f)

end

function build_chains(ps,Cs)

    # dimensions of the new C matrix
    nbf5_total = 0
    ndoc_total = 0
    nsoc_total = 0
    nbf_total = 0
    for p in ps
        nbf5_total += p.nbf5
        ndoc_total += p.ndoc
        nsoc_total += p.nsoc
        nbf_total += p.nbf
    end
    ndns_total = ndoc_total + nsoc_total

    # new C matrix
    C_new = zeros(nbf_total,nbf_total)

    # indices of the new C matrix
    idoc_new = 1
    icwo_new = nbf5_total
    isoc_new = ndoc_total+1
    iout = nbf5_total+1
    ibf = 1
    for (p,C) in zip(ps,Cs)
        
        # add double occupied (and coupled) orbitals
        for i in 1:p.ndoc	
            idx = p.no1 + i

	    C_new[ibf:ibf+p.nbf-1,idoc_new] = C[:,idx]
	    idoc_new += 1

            ll = p.no1 + p.ndns + p.ncwo*(p.ndoc - i) + 1
            ul = p.no1 + p.ndns + p.ncwo*(p.ndoc - i + 1)
            C_new[ibf:ibf+p.nbf-1,icwo_new-(ul-ll):icwo_new] = C[:,ll:ul]

	    icwo_new += -(ul-ll + 1)
        end

	# add single occupied orbitals
        C_new[ibf:ibf+p.nbf-1,isoc_new:isoc_new+p.nsoc-1] = C[:,p.ndoc+1:p.ndns]
	isoc_new += p.nsoc

	# add unoccupied orbitals 
	C_new[ibf:ibf+p.nbf-1,iout:iout+p.nbf-p.nbf5-1] = C[:,p.nbf5+1:p.nbf]
	iout += p.nbf-p.nbf5

	ibf += p.nbf
    end

    return C_new

end

function guess_gamma_trigonometric(p)
    gamma = zeros(p.nv)
    for i in 1:p.ndoc
        gamma[i] = acos(sqrt(2.0*0.999-1.0))
        for j in 1:p.ncwo-1
            ig = p.ndoc+(i-1)*(p.ncwo-1)+j
            gamma[ig] = asin(sqrt(1.0/(p.ncwo-j+1)))
        end
    end
    return gamma
end
