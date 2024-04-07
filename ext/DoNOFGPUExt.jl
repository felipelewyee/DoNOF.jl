module DoNOFGPUExt

using DoNOF

using CUDA, KernelAbstractions, NNlibCUDA, cuTENSOR

using LoopVectorization
using Tullio
using TensorOperations

using LinearAlgebra
using Printf

using Optim
using SparseArrays
using GaussianBasis
using JLD

USE_CUDA = true

#function DoNOF.compute_integrals(bset,p)
#
#    # Overlap, Kinetics, Potential
#    S = overlap(bset)
#    T = kinetic(bset)
#    V = nuclear(bset)
#    H = T + V
#    I = Float64[]
#    b_mnl = Float64[]
#    if (!p.RI)
#        # Integrales de Repulsión Electrónica, ERIs (mu nu | sigma lambda)
#	I = ERI_2e4c(bset)
#    else
#        if p.spherical
#	    aux = BasisSet(bset.name*"-jkfit",bset.atoms)
#        else
#	    aux = BasisSet(bset.name*"-jkfit",bset.atoms,spherical=false,lib=:acsint)
#        end
#        Iaux = ERI_2e3c(bset,aux)
#        G = ERI_2e2c(aux)
#
#        evals,evecs = eigen(G)
#        sqrtinv = Float64[]
#        for i in 1:size(evals)[1]
#            if (evals[i]<0.0)
#                append!(sqrtinv,0.0)
#            else
#                append!(sqrtinv,1/sqrt(evals[i]))
#            end
#        end
#        Gmsqrt = evecs*Diagonal(sqrtinv)*evecs'
#        @tullio b_mnl[m,n,l] := Iaux[m,n,k]*Gmsqrt[k,l]
#
#        p.nbfaux = size(b_mnl)[3]
#    end
#
#    if (!p.RI)
#        I = CuArray(I)
#    else
#       if(p.gpu_bits==64)
#            b_mnl = CuArray(b_mnl)
#       else
#            b_mnl = cu(b_mnl)
#       end
#    end
#    println("test")
#    return S,T,V,H,I,b_mnl
#
#end

end
