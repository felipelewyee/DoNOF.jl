module DoNOF

using LoopVectorization
using Tullio
using TensorOperations

using LinearAlgebra
using Printf
using SpecialFunctions

using Optim
using GaussianBasis
using FileIO, JLD2

using FiniteDiff

include("io.jl")
include("param.jl")
include("integrals.jl")
include("utils.jl")
include("pnof.jl")
include("minimization.jl")
include("guess.jl")
include("postpnof.jl")
include("energy.jl")
include("output.jl")
include("grads.jl")
include("hess.jl")
include("optgeo.jl")

using PrecompileTools

@setup_workload begin
mol = """
0 1
 O  0.0000   0.000   0.121
 H  0.0000   0.751  -0.485
 H  0.0000  -0.751  -0.485
"""

@compile_workload begin
bset,p = DoNOF.molecule(mol,"cc-pvdz",spherical=true)

p.ipnof = 8
p.ista = 0

p.maxit = 2

DoNOF.energy(bset,p)
end

end

end # module
