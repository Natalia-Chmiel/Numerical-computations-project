using LinearAlgebra

function solve_FDM(V::Vector{Float64}, x::Vector{Float64}, m::Float64, planck::Float64)
    N = length(x)
    dx = x[2] - x[1]
    
    H = zeros(N, N)
    
    for i in 1:N
        H[i, i] = V[i] + planck^2 / (m * dx^2) 
        if i > 1
            H[i, i-1] = (-1) * planck^2 / (2 * m * dx^2) 
        end
        if i < N
            H[i, i+1] = (-1) * planck^2 / (2 * m * dx^2)
        end
    end
    
    eigenvalues, eigenvectors = eigen(H)
    
    return eigenvalues, eigenvectors
end