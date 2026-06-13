using LinearAlgebra


"""
solve_FDM(V::Vector{Float64}, x::Vector{Float64}, m::Float64, planck::Float64)

Rozwiązuje jednowymiarowe, niezależne od czasu równanie Schrodingera przy użyciu Metody Różnic Skończonych (FDM).

# Argumenty:
- `V`: Wektor (`Vector{Float64}`) zawierający wartości potencjału w węzłach siatki przestrzennej.
- `x`: Wektor (`Vector{Float64}`) punktów siatki przestrzennej.
- `m`: Masa cząstki (`Float64`).
- `planck`: Zredukowana stała Plancka (`Float64`).

# Zwraca:
- Krotkę `(eigenvalues, eigenvectors)`:
  - `eigenvalues`: Wektor obliczonych wartości własnych (poziomów energii).
  - `eigenvectors`: Macierz wektorów własnych (funkcji falowych), gdzie każda kolumna odpowiada danej energii. Wektory są automatycznie znormalizowane.
"""
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