function compute_gamma(p)
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

function read_C(;title = "donof")

    C = load(title*".jld")["C"]
    return C

end

function read_gamma(;title = "donof")

    gamma = load(title*".jld")["gamma"]
    return gamma

end

function read_fmiug0(;title = "donof")
   
    fmiug0 = load(title*".jld")["fmiug0"]
    return fmiug0

end

function read_all(;title = "donof")
    C = read_C(title=title)
    gamma = read_gamma(title=title)
    fmiug0 = read_fmiug0(title=title)

    return C,gamma,fmiug0

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
