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

    gamma = np.load(title*"jld")["gamma"]
    return gamma

end

function read_fmiug0(;title = "donof")
   
    fmiug0 = load(title*"jld")["fmiug0"]
    return fmiug0

end

function read_all(;title = "donof")
    C = read_C(title=title)
    gamma = read_gamma(title=title)
    fmiug0 = read_fmiug0(title=title)

    return C,gamma,fmiug0

end

