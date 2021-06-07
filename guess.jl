module guess

function compute_gamma(p)
    gamma = zeros(p.nbf5)
    for i in 1:p.ndoc
        gamma[i] = acos(sqrt(2.0*0.999-1.0))
        for j in 1:p.ncwo-1
		ig = p.ndoc+(i-1)*(p.ncwo-1)+j
	    gamma[ig] = asin(sqrt(1.0/(p.ncwo-j+1)))

        end
    end
    return gamma
end

end
