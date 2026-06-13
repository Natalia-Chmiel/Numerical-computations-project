using LinearAlgebra

function solve_FDM(V::Vector{Float64}, x::Vector{Float64}, m::Float64, planck::Float64)

    N = length(x)
    dx = x[2] - x[1]

    C = planck^2 / (2.0 * m * dx^2) 

    main_diag = V .+ (2.0 * C) 
    under_diag = fill(-C, N - 1)

    H = SymTridiagonal(main_diag, under_diag)
    
    eigenvalues, eigenvectors = eigen(H)
    eigenvectors = eigenvectors ./ sqrt.(dx) 
    
    return eigenvalues, eigenvectors
end