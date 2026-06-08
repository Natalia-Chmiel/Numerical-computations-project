module SchrodingerSolver

using LinearAlgebra

include("fdm.jl")
include("numerov.jl")

export solve_FDM, solve_Numerov

end