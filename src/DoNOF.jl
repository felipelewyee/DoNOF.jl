module DoNOF

using LoopVectorization
using CUDA, CUDAKernels, KernelAbstractions
using Tullio

using LinearAlgebra
using Printf

using Optim
using SparseArrays
using GaussianBasis
using JLD

include("io.jl")
include("param.jl")
include("integrals.jl")
include("utils.jl")
include("pnof.jl")
include("minimization.jl")
include("guess.jl")
include("postpnof.jl")
include("energy.jl")

mol = """
  O    0.0000000    0.0099701   -0.0000000
  H   -0.0000000   -0.3962400    0.7727381
  H    0.0000000   -0.3962400   -0.7727381
"""

bset,p = DoNOF.molecule(mol,"cc-pvdz",spherical=false)

p.ipnof = 8

p.RI = true
p.gpu = false

p.method = "ID"

p.maxitid = 2
p.maxit = 2

E,C,gamma,fmiug0 = DoNOF.energy(bset,p,do_hfidr=true,do_m_diagnostic=true)

end # module
